import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

# --- 1. CONFIGURATION ---
st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "type": "Circular", "expected_genes": 4300},
    "Staphylococcus aureus": {"ref_gc": 32.8, "type": "Circular", "expected_genes": 2800},
    "Homo sapiens (Partial mRNA)": {"ref_gc": 41.0, "type": "Linear", "expected_genes": 1}
}

# --- 2. BIOLOGICAL FUNCTIONS ---
def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def remove_adapters(reads, adapter_seq, min_keep_len):
    cleaned_reads = []
    for read in reads:
        if adapter_seq and adapter_seq in read:
            read = read.split(adapter_seq)[0]
        if len(read) >= min_keep_len:
            cleaned_reads.append(read)
    return cleaned_reads

def calculate_n50(lengths):
    lengths.sort(reverse=True)
    total_sum = sum(lengths)
    current_sum = 0
    for length in lengths:
        current_sum += length
        if current_sum >= total_sum / 2:
            return length
    return 0

def find_all_orfs(sequence, min_len=300, allow_partial=True):
    found_genes = []
    pattern_str = r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA%s))' % (min_len // 3, '|$' if allow_partial else '')
    pattern = re.compile(pattern_str)
    
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                if len(gene_seq) % 3 != 0:
                    gene_seq = gene_seq[:-(len(gene_seq) % 3)]
                if len(gene_seq) < min_len: continue
                
                start_pos = match.start() + frame
                found_genes.append({
                    "Strand": strand, "Start": int(start_pos), "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)), "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2),
                    "Sequence": gene_seq
                })
    
    sorted_genes = sorted(found_genes, key=lambda x: x['Length'], reverse=True)
    final_genes, covered = [], []
    for g in sorted_genes:
        if not any(max(g['Start'], s) < min(g['End'], e) for s, e in covered):
            final_genes.append(g)
            covered.append((g['Start'], g['End']))
    
    final_genes = sorted(final_genes, key=lambda x: x['Start'])
    for i, gene in enumerate(final_genes):
        gene['Name'] = f"ORF_{i+1}"
        
    return final_genes

# --- 3. UI SIDEBAR ---
st.sidebar.header("⚙️ Pipeline Settings")
viz_mode = st.sidebar.radio("Map Visualization:", ("Linear Track", "Circular Map"))
min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 50, 1000, 300)
allow_partial = st.sidebar.checkbox("Allow Partial Genes", value=True)

st.sidebar.markdown("---")
st.sidebar.subheader("🛡️ QC & Trimming")
adapter_seq = st.sidebar.text_input("Adapter Sequence", "AGATCGGAAGAG")
min_read_len = st.sidebar.slider("Min Length Threshold", 10, 500, 30)

# --- 4. PROCESSING ---
st.title("🧬 De Novo: Professional Genomic Suite")

uploaded_file = st.file_uploader("Upload FASTA or FASTQ", type=["fasta", "fa", "fastq", "fq", "gz", "txt"])

if uploaded_file:
    try:
        content = (gzip.decompress(uploaded_file.read()).decode("utf-8") 
                   if uploaded_file.name.endswith('.gz') else uploaded_file.read().decode("utf-8"))
        lines = content.splitlines()
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        # Raw Data Metrics
        if is_fasta:
            raw_reads = []
            curr = []
            for l in lines:
                if l.startswith(">"):
                    if curr: raw_reads.append("".join(curr))
                    curr = []
                else: curr.append(l.strip())
            if curr: raw_reads.append("".join(curr))
        else:
            raw_reads = [l.strip() for l in lines[1::4] if l.strip()]

        raw_count = len(raw_reads)
        raw_lengths = [len(r) for r in raw_reads]
        
        st.info(f"**Format:** {'FASTA' if is_fasta else 'FASTQ'} | **Initial Reads:** {raw_count}")

        if st.button("🚀 Run Full Pipeline"):
            # Trimming Logic
            processed_reads = remove_adapters(raw_reads, adapter_seq, min_read_len)
            proc_count = len(processed_reads)
            proc_lengths = [len(r) for r in processed_reads]
            
            if proc_count == 0:
                st.error("QC filtered out all reads. Check your length threshold.")
                st.stop()

            full_seq = "NNNNN".join(processed_reads)
            total_len = len(full_seq)
            sample_gc = round((full_seq.count('G') + full_seq.count('C')) / (total_len - full_seq.count('N') + 1) * 100, 2)
            
            raw_genes = find_all_orfs(full_seq, min_len=min_orf_len, allow_partial=allow_partial)
            df = pd.DataFrame(raw_genes)
            
            t1, t2, t3 = st.tabs(["📊 QC Report", "🏗️ Assembly Metrics", "🧬 Genomic Map"])
            
            with t1:
                st.subheader("🛡️ Trimming & QC Comparison")
                qc_df = pd.DataFrame({
                    "Metric": ["Total Reads/Contigs", "Average Length (bp)", "Max Length (bp)"],
                    "Before (Raw)": [raw_count, f"{np.mean(raw_lengths):.1f}", max(raw_lengths)],
                    "After (Processed)": [proc_count, f"{np.mean(proc_lengths):.1f}", max(proc_lengths)],
                    "Change": [proc_count - raw_count, f"{np.mean(proc_lengths)-np.mean(raw_lengths):.1f}", max(proc_lengths)-max(raw_lengths)]
                })
                st.table(qc_df)

            with t2:
                st.subheader("🏗️ Assembly Quality Statistics")
                c1, c2, c3 = st.columns(3)
                
                # Calculate N50
                n50_val = calculate_n50(proc_lengths)
                coding_density = (df['Length'].sum() / total_len) * 100 if not df.empty else 0
                
                c1.metric("N50 Score", f"{n50_val} bp")
                c2.metric("Coding Density", f"{coding_density:.1f}%")
                c3.metric("GC Content", f"{sample_gc}%")
                
                if not df.empty:
                    st.write("#### ORF Feature Table")
                    st.dataframe(df.drop(columns=['Sequence']), use_container_width=True)

            with t3:
                if df.empty:
                    st.warning("No genes found.")
                elif viz_mode == "Linear Track":
                    fig = go.Figure()
                    for _, row in df.iterrows():
                        clr = "#00CC96" if row['Strand'] == "Forward" else "#EF553B"
                        fig.add_trace(go.Bar(name=row['Name'], x=[row['Length']], y=[row['Strand']], base=[row['Start']], orientation='h', marker_color=clr))
                    fig.update_layout(template="plotly_dark", barmode='stack', title="Linear ORF Map")
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    df['S_Ang'], df['E_Ang'] = (df['Start']/total_len)*360, (df['End']/total_len)*360
                    fig = go.Figure()
                    for _, row in df.iterrows():
                        track, clr = (2.1, "#00CC96") if row['Strand']=="Forward" else (1.6, "#EF553B")
                        fig.add_trace(go.Barpolar(name=row['Name'], r=[0.4], theta=[(row['S_Ang']+row['E_Ang'])/2], width=[max(1, row['E_Ang']-row['S_Ang'])], base=track, marker_color=clr))
                    fig.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(visible=False)), height=700, title="Circular Genome Map")
                    st.plotly_chart(fig, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

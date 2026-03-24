import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

# --- 1. CONFIGURATION ---
st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "type": "Circular"},
    "Staphylococcus aureus": {"ref_gc": 32.8, "type": "Circular"},
    "Homo sapiens (Partial mRNA)": {"ref_gc": 41.0, "type": "Linear"}
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
st.sidebar.subheader("🛡️ Trimming Settings")
trim_active = st.sidebar.toggle("Enable QC Trimming", value=True)
adapter_seq = st.sidebar.text_input("Adapter to Remove", "AGATCGGAAGAG")
min_read_len = st.sidebar.slider("Min Read Length to Keep", 10, 100, 30)

# --- 4. PROCESSING ---
st.title("🧬 De Novo: Professional Genomic Suite")

uploaded_file = st.file_uploader("Upload FASTA or FASTQ", type=["fasta", "fa", "fastq", "fq", "gz", "txt"])

if uploaded_file:
    try:
        content = (gzip.decompress(uploaded_file.read()).decode("utf-8") 
                   if uploaded_file.name.endswith('.gz') else uploaded_file.read().decode("utf-8"))
        lines = content.splitlines()
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        # Raw Data Extraction
        if is_fasta:
            raw_reads = ["".join([l.strip() for l in lines if l.strip() and not l.startswith('>')])]
        else:
            raw_reads = [l.strip() for l in lines[1::4] if l.strip()]

        # 1. Capture Raw Metrics
        raw_count = len(raw_reads)
        raw_avg_len = sum(len(r) for r in raw_reads) / raw_count if raw_count > 0 else 0

        if st.button("🚀 Run Analysis"):
            # 2. Apply Trimming/QC
            if trim_active:
                processed_reads = remove_adapters(raw_reads, adapter_seq, min_read_len)
            else:
                processed_reads = raw_reads

            # 3. Capture Processed Metrics
            proc_count = len(processed_reads)
            proc_avg_len = sum(len(r) for r in processed_reads) / proc_count if proc_count > 0 else 0
            
            if proc_count == 0:
                st.error("QC filters removed all reads. Please adjust settings.")
                st.stop()

            full_seq = "NNNNN".join(processed_reads[:1000])
            total_len = len(full_seq)
            sample_gc = round((full_seq.count('G') + full_seq.count('C')) / (total_len - full_seq.count('N') + 1) * 100, 2)
            
            raw_genes = find_all_orfs(full_seq, min_len=min_orf_len, allow_partial=allow_partial)
            df = pd.DataFrame(raw_genes)
            
            t1, t2 = st.tabs(["📊 QC & Metrics", "🧬 Genomic Map"])
            
            with t1:
                st.subheader("🛡️ Quality Control Comparison")
                qc_col1, qc_col2, qc_col3 = st.columns(3)
                qc_col1.metric("Reads: Before vs After", f"{proc_count}", f"{proc_count - raw_count}")
                qc_col2.metric("Avg Length: Before vs After", f"{proc_avg_len:.1f}bp", f"{proc_avg_len - raw_avg_len:.1f}bp")
                qc_col3.metric("Retention Rate", f"{(proc_count/raw_count)*100:.1f}%")

                st.markdown("---")
                st.write("#### Detailed Comparison Table")
                comparison_df = pd.DataFrame({
                    "Stage": ["Raw Data (Input)", "Processed Data (Output)"],
                    "Total Reads": [raw_count, proc_count],
                    "Avg Read Length": [f"{raw_avg_len:.1f} bp", f"{proc_avg_len:.1f} bp"],
                    "Max Read Length": [max(len(r) for r in raw_reads), max(len(r) for r in processed_reads)]
                })
                st.table(comparison_df)

                st.markdown("---")
                st.metric("Final ORFs Detected", len(df))
                if not df.empty:
                    st.dataframe(df.drop(columns=['Sequence']), use_container_width=True)

            with t2:
                if df.empty:
                    st.warning("No genes found.")
                elif viz_mode == "Linear Track":
                    fig = go.Figure()
                    for _, row in df.iterrows():
                        clr = "#00CC96" if row['Strand'] == "Forward" else "#EF553B"
                        fig.add_trace(go.Bar(
                            name=row['Name'], x=[row['Length']], y=[row['Strand']], 
                            base=[row['Start']], orientation='h', marker_color=clr,
                            hovertext=f"Name: {row['Name']}<br>Start: {row['Start']}<br>GC: {row['GC %']}%",
                            hoverinfo="text"
                        ))
                    fig.update_layout(template="plotly_dark", barmode='stack', title="Linear ORF Map")
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    df['S_Ang'], df['E_Ang'] = (df['Start']/total_len)*360, (df['End']/total_len)*360
                    fig = go.Figure()
                    for _, row in df.iterrows():
                        track, clr = (2.1, "#00CC96") if row['Strand']=="Forward" else (1.6, "#EF553B")
                        fig.add_trace(go.Barpolar(
                            name=row['Name'], r=[0.4], theta=[(row['S_Ang']+row['E_Ang'])/2], 
                            width=[max(1, row['E_Ang']-row['S_Ang'])], base=track, marker_color=clr,
                            hovertext=f"Name: {row['Name']}<br>Start: {row['Start']}<br>GC: {row['GC %']}%",
                            hoverinfo="text"
                        ))
                    fig.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(visible=False)), height=700, title="Circular Genome Map")
                    st.plotly_chart(fig, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

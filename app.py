import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

# --- 1. CONFIGURATION ---
st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200},
    "Homo sapiens (Fragment)": {"ref_gc": 41.0, "expected_genes": 20000}
}

# --- 2. BIOLOGICAL FUNCTIONS ---
def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=300, allow_partial=True):
    found_genes = []
    # Regex: Start with ATG, triplets, end with Stop OR End of string ($)
    pattern_str = r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA%s))' % (min_len // 3, '|$' if allow_partial else '')
    pattern = re.compile(pattern_str)
    
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                # Clean tailing bases
                if len(gene_seq) % 3 != 0:
                    gene_seq = gene_seq[:-(len(gene_seq) % 3)]
                
                if len(gene_seq) < min_len: continue
                
                start_pos = match.start() + frame
                found_genes.append({
                    "Strand": strand, 
                    "Start": int(start_pos), 
                    "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)), 
                    "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2),
                    "Sequence": gene_seq
                })
    
    if not found_genes: return []
    
    # FILTER: Remove excessive overlaps to fix "149% Coding Density"
    # We sort by length and keep the longest ORFs that don't overlap significantly
    sorted_genes = sorted(found_genes, key=lambda x: x['Length'], reverse=True)
    final_genes = []
    covered_regions = [] # List of (start, end)

    for gene in sorted_genes:
        is_overlapping = False
        for start, end in covered_regions:
            # If more than 50% of the gene overlaps an existing one, discard it
            overlap = max(0, min(gene['End'], end) - max(gene['Start'], start))
            if overlap > (gene['Length'] * 0.5):
                is_overlapping = True
                break
        if not is_overlapping:
            final_genes.append(gene)
            covered_regions.append((gene['Start'], gene['End']))
            
    return final_genes

# --- 3. UI SIDEBAR ---
st.sidebar.header("⚙️ Pipeline Settings")
min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 50, 2000, 300)
allow_partial = st.sidebar.checkbox("Allow Partial Genes", value=True)
read_limit = st.sidebar.slider("Read Limit (FASTQ only)", 100, 5000, 1000)

# --- 4. DATA PROCESSING ---
st.title("🧬 De Novo: Professional Genomic Suite")

uploaded_file = st.file_uploader("Upload FASTA or FASTQ", type=["fasta", "fa", "fastq", "fq", "gz", "txt"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.gz'):
            content = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            content = uploaded_file.read().decode("utf-8")
        
        lines = content.splitlines()
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        if is_fasta:
            # Stitched FASTA logic
            full_seq = "".join([l.strip() for l in lines if l.strip() and not l.startswith('>')])
            raw_data = [full_seq]
        else:
            # FASTQ logic
            raw_data = [l.strip() for l in lines[1::4] if l.strip()]

        final_scaffold = "NNNNN".join(raw_data[:read_limit])
        total_len = len(final_scaffold)
        sample_gc = round((final_scaffold.count('G') + final_scaffold.count('C')) / (total_len - final_scaffold.count('N')) * 100, 2)
        closest = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))

        st.info(f"**Format:** {'FASTA' if is_fasta else 'FASTQ'} | **Scaffold Length:** {total_len} bp")
        st.success(f"**Predicted Species:** {closest} (GC: {sample_gc}%)")

        if st.button("🚀 Run Full Analysis"):
            genes = find_all_orfs(final_scaffold, min_len=min_orf_len, allow_partial=allow_partial)
            df = pd.DataFrame(genes)
            
            tab1, tab2, tab3 = st.tabs(["📊 Metrics", "🏗️ GC Skew", "🧬 Circular Map"])
            
            with tab1:
                col1, col2 = st.columns(2)
                col1.metric("ORFs Detected", len(df))
                coding_len = df['Length'].sum() if not df.empty else 0
                col2.metric("Coding Density", f"{round((coding_len/total_len)*100, 1)}%")
                
                if not df.empty:
                    st.subheader("🔍 Top ORF for BLAST")
                    st.text_area("Copy/Paste into NCBI:", df.sort_values('Length', ascending=False).iloc[0]['Sequence'], height=100)
                    st.dataframe(df.drop(columns=['Sequence']), use_container_width=True)

            with tab2:
                # GC Skew Logic
                win = 500
                p_skew, skews = [], []
                clean_dna = final_scaffold.replace("N", "")
                for i in range(0, len(clean_dna) - win, win):
                    sub = clean_dna[i:i+win]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g-c)/(g+c) if (g+c)>0 else 0)
                    p_skew.append(i)
                fig_s = go.Figure()
                s_np = np.array(skews)
                fig_s.add_trace(go.Scatter(x=p_skew, y=np.where(s_np>=0, s_np, 0), fill='tozeroy', line_color='#00CC96', name='G>C'))
                fig_s.add_trace(go.Scatter(x=p_skew, y=np.where(s_np<0, s_np, 0), fill='tozeroy', line_color='#EF553B', name='C>G'))
                fig_s.update_layout(template="plotly_dark", height=400, title="GC Skew")
                st.plotly_chart(fig_s, use_container_width=True)

            with tab3:
                if not df.empty:
                    df['S_Ang'] = (df['Start'] / total_len) * 360
                    df['E_Ang'] = (df['End'] / total_len) * 360
                    fig_c = go.Figure()
                    for _, row in df.iterrows():
                        track = 2.1 if row['Strand'] == "Forward" else 1.6
                        clr = "#00CC96" if row['Strand'] == "Forward" else "#EF553B"
                        fig_c.add_trace(go.Barpolar(r=[0.4], theta=[(row['S_Ang']+row['E_Ang'])/2], 
                                                    width=[max(1, row['E_Ang']-row['S_Ang'])], 
                                                    base=track, marker_color=clr))
                    fig_c.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(visible=False)), height=600)
                    st.plotly_chart(fig_c, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

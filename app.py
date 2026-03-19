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
    "Homo sapiens (Partial Fragment)": {"ref_gc": 41.0, "expected_genes": 20000}
}

# --- 2. BIOLOGICAL ENGINE ---
def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=150, allow_partial=True):
    found_genes = []
    # Logic: Start with ATG, follow triplets, end with Stop OR the end of the string ($)
    if allow_partial:
        pattern = re.compile(r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA|$))' % (min_len // 3))
    else:
        pattern = re.compile(r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA))' % (min_len // 3))
    
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                # Ensure the sequence length is a multiple of 3 (triplets)
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
    return found_genes

# --- 3. SIDEBAR ---
st.sidebar.header("⚙️ Analysis Settings")
min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 50, 1000, 150)
allow_partial = st.sidebar.checkbox("Allow Partial Genes (No Stop Codon)", value=True)

# --- 4. UPLOAD & PROCESSING ---
st.title("🧬 De Novo: Professional Genomic Suite")
uploaded_file = st.file_uploader("Upload FASTA or FASTQ", type=["fasta", "fa", "fastq", "fq", "txt"])

if uploaded_file:
    try:
        content = uploaded_file.read().decode("utf-8")
        lines = content.splitlines()
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        if is_fasta:
            # STITCHING: Joins all sequence lines into one continuous string
            full_seq = "".join([l.strip() for l in lines if l.strip() and not l.startswith('>')])
            raw_data = [full_seq]
        else:
            raw_data = [l.strip() for l in lines[1::4] if l.strip()]

        total_dna = "".join(raw_data)
        sample_gc = round((total_dna.count('G') + total_dna.count('C')) / len(total_dna) * 100, 2)
        closest = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))

        st.info(f"**Format:** {'FASTA' if is_fasta else 'FASTQ'} | **Length:** {len(total_dna)} bp | **Predicted:** {closest}")

        if st.button("🚀 Run Genomic Analysis"):
            genes = find_all_orfs(total_dna, min_len=min_orf_len, allow_partial=allow_partial)
            df = pd.DataFrame(genes)
            
            tab1, tab2 = st.tabs(["📊 Data Metrics", "🧬 Circular Annotation"])
            
            with tab1:
                st.metric("Total ORFs Found", len(df))
                if not df.empty:
                    st.subheader("🔍 Top Sequence for BLAST")
                    top_seq = df.sort_values('Length', ascending=False).iloc[0]['Sequence']
                    st.text_area("Copy/Paste into NCBI BLAST:", value=top_seq, height=150)
                    st.dataframe(df.drop(columns=['Sequence']), use_container_width=True)

            with tab2:
                if not df.empty:
                    df['S_Ang'] = (df['Start'] / len(total_dna)) * 360
                    df['E_Ang'] = (df['End'] / len(total_dna)) * 360
                    
                    fig = go.Figure()
                    for _, row in df.iterrows():
                        track = 2.2 if row['Strand'] == "Forward" else 1.7
                        clr = "#00CC96" if row['Strand'] == "Forward" else "#EF553B"
                        fig.add_trace(go.Barpolar(
                            r=[0.4], theta=[(row['S_Ang']+row['E_Ang'])/2],
                            width=[max(2, row['E_Ang']-row['S_Ang'])], 
                            base=track, marker_color=clr, name=row['Strand']
                        ))
                    fig.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(visible=False)), height=600)
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("No genes found. Try lowering the Minimum ORF Length.")

    except Exception as e:
        st.error(f"Error: {e}")

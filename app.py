import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

# --- 1. CONFIGURATION ---
st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")

# --- 2. BIOLOGICAL FUNCTIONS ---
def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

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
    
    # Filter Overlaps (Keeping only the most significant ORF per locus to avoid clutter)
    sorted_genes = sorted(found_genes, key=lambda x: x['Length'], reverse=True)
    final_genes, covered = [], []
    for g in sorted_genes:
        if not any(max(g['Start'], s) < min(g['End'], e) for s, e in covered):
            final_genes.append(g)
            covered.append((g['Start'], g['End']))
    
    # Return sorted by original genomic position
    return sorted(final_genes, key=lambda x: x['Start'])

# --- 3. UI SIDEBAR ---
st.sidebar.header("⚙️ Pipeline Settings")
min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 50, 1000, 300)
allow_partial = st.sidebar.checkbox("Allow Partial Genes", value=True)
# Added toggle for Gapless mode
map_style = st.sidebar.radio("Map Style:", ("Gapless (Stitched)", "Real-Coordinates (With Gaps)"))

# --- 4. PROCESSING ---
st.title("🧬 De Novo: Professional Genomic Suite")

uploaded_file = st.file_uploader("Upload FASTA or FASTQ", type=["fasta", "fa", "fastq", "fq", "gz", "txt"])

if uploaded_file:
    try:
        content = (gzip.decompress(uploaded_file.read()).decode("utf-8") 
                   if uploaded_file.name.endswith('.gz') else uploaded_file.read().decode("utf-8"))
        lines = content.splitlines()
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        full_seq = "".join([l.strip() for l in lines if l.strip() and not l.startswith('>')]) if is_fasta else "NNNNN".join([l.strip() for l in lines[1::4] if l.strip()][:1000])
        total_len = len(full_seq)

        if st.button("🚀 Run Analysis"):
            genes = find_all_orfs(full_seq, min_len=min_orf_len, allow_partial=allow_partial)
            df = pd.DataFrame(genes)
            
            t1, t2 = st.tabs(["📊 Data Metrics", "🧬 Linear Map"])
            
            with t1:
                st.metric("Total Genes Found", len(df))
                if not df.empty:
                    st.dataframe(df.drop(columns=['Sequence']), use_container_width=True)

            with t2:
                if df.empty:
                    st.warning("No genes found. Try decreasing the Minimum ORF Length.")
                else:
                    fig = go.Figure()
                    for strand in ["Forward", "Reverse"]:
                        sdf = df[df["Strand"] == strand].copy()
                        clr = "#00CC96" if strand == "Forward" else "#EF553B"
                        
                        if map_style == "Gapless (Stitched)":
                            # Remove gaps by calculating cumulative length per strand
                            lengths = sdf["Length"].tolist()
                            # Cumulative sum starting at 0
                            cumulative_bases = [0] + list(np.cumsum(lengths))[:-1]
                            
                            fig.add_trace(go.Bar(
                                x=lengths, 
                                y=[strand]*len(sdf), 
                                base=cumulative_bases, 
                                orientation='h', 
                                marker_color=clr, 
                                name=strand,
                                hovertext=[f"Original Position: {s}bp" for s in sdf["Start"]]
                            ))
                        else:
                            # Original behavior with genomic gaps
                            fig.add_trace(go.Bar(
                                x=sdf["Length"], 
                                y=[strand]*len(sdf), 
                                base=sdf["Start"], 
                                orientation='h', 
                                marker_color=clr, 
                                name=strand
                            ))

                    fig.update_layout(
                        template="plotly_dark", 
                        barmode='overlay', # Overlay allows Forward/Reverse to sit tightly
                        height=400, 
                        xaxis_title="Coding Position (bp)" if map_style == "Gapless (Stitched)" else "Genomic Position (bp)",
                        yaxis=dict(categoryorder='array', categoryarray=['Forward', 'Reverse'])
                    )
                    st.plotly_chart(fig, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

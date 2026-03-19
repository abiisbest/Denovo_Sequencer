import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

# --- 1. CONFIGURATION & DATABASE ---
st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200},
    "Homo sapiens (Partial Fragment)": {"ref_gc": 41.0, "expected_genes": 20000}
}

# --- 2. CORE BIOLOGICAL FUNCTIONS ---
def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=300, allow_partial=True):
    found_genes = []
    # Regex: ATG, triplets, then Stop OR End of string ($)
    if allow_partial:
        pattern = re.compile(r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA|$))' % (min_len // 3))
    else:
        pattern = re.compile(r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA))' % (min_len // 3))
    
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                # Clean up tailing bases that don't form a full triplet
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

# --- 3. SIDEBAR SETTINGS ---
st.sidebar.header("⚙️ Pipeline Parameters")
read_limit = st.sidebar.slider("Max Reads/Sequences to Process", 1, 5000, 1000)
min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 50, 2000, 150)
allow_partial = st.sidebar.checkbox("Allow Partial Genes (No Stop Codon)", value=True)

st.sidebar.subheader("🛡️ QC & Trimming")
enable_trim = st.sidebar.checkbox("Enable End Trimming", value=False)
trim_val = st.sidebar.slider("Trim from ends (bp)", 0, 50, 5)

# --- 4. DATA UPLOAD & INITIAL PROCESSING ---
st.title("🧬 De Novo: Professional Genomic Suite")
st.markdown("---")

uploaded_file = st.file_uploader("Upload FASTA, FASTQ, or GZ", type=["fasta", "fa", "fastq", "fq", "gz", "txt"])

if uploaded_file:
    try:
        # Handle Compression
        if uploaded_file.name.endswith('.gz'):
            content = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            content = uploaded_file.read().decode("utf-8")
        
        lines = content.splitlines()
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        if is_fasta:
            # STITCHING: Combine all non-header lines into one scaffold
            full_seq = "".join([l.strip() for l in lines if l.strip() and not l.startswith('>')])
            raw_data = [full_seq]
            format_label = "FASTA (Stitched)"
        else:
            # FASTQ: Every 2nd line in 4-line blocks
            raw_data = [l.strip() for l in lines[1::4] if l.strip()]
            format_label = "FASTQ (Multi-read)"

        # Species Prediction based on GC
        total_dna_string = "".join(raw_data[:100])
        sample_gc = round((total_dna_string.count('G') + total_dna_string.count('C')) / len(total_dna_string) * 100, 2)
        closest = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))
        
        st.subheader("🤖 Identification Results")
        st.info(f"**Format:** {format_label} | **Sample GC:** {sample_gc}%")
        st.success(f"**Predicted Species:** {closest}")
        
        # --- 5. EXECUTION ENGINE ---
        if st.button("🚀 Run Analysis Pipeline"):
            # Step 1: Cleaning/Filtering
            if enable_trim:
                processed = [s[trim_val:-trim_val] for s in raw_data if len(s) > (min_orf_len + 2*trim_val)]
            else:
                processed = [s for s in raw_data if len(s) >= (min_orf_len // 2)]
            
            # Step 2: Assembly (Pseudo-Scaffold)
            final_scaffold = "NNNNN".join(processed[:read_limit])
            total_len = len(final_scaffold)
            
            # Step 3: Annotation
            orfs = find_all_orfs(final_scaffold, min_len=min_orf_len, allow_partial=allow_partial)
            genes_df = pd.DataFrame(orfs)
            if not genes_df.empty:
                genes_df = genes_df.sort_values('Start').drop_duplicates(subset=['Start'])

            # --- 6. VISUALIZATION TABS ---
            tab1, tab2, tab3 = st.tabs(["📊 QC & Metrics", "🏗️ Structural Analysis", "🧬 Circular Map"])

            with tab1:
                st.subheader("🛡️ Read Quality & Yield")
                col_m1, col_m2 = st.columns(2)
                col_m1.metric("Input Fragments", len(raw_data))
                col_m2.metric("Valid Sequences", len(processed), f"{len(processed)-len(raw_data)}")
                
                results_summary = {
                    "Metric": ["Scaffold Length", "GC Content", "ORFs Detected", "Coding Density"],
                    "Value": [f"{total_len} bp", f"{sample_gc}%", f"{len(genes_df)}", f"{round((genes_df['Length'].sum()/total_len)*100,1) if not genes_df.empty else 0}%"]
                }
                st.table(pd.DataFrame(results_summary))

            with tab2:
                st.subheader("🏗️ GC Skew Analysis (Genomic Landmarks)")
                window = 500
                # Filter out Ns for a cleaner skew calculation
                clean_seq = final_scaffold.replace("N", "")
                p_skew, skews = [], []
                for i in range(0, len(clean_seq) - window, window):
                    sub = clean_seq[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)

                fig_skew = go.Figure()
                skews_np = np.array(skews)
                fig_skew.add_trace(go.Scatter(x=p_skew, y=np.where(skews_np >= 0, skews_np, 0), fill='tozeroy', line_color='#00CC96', name='Positive Skew (G > C)'))
                fig_skew.add_trace(go.Scatter(x=p_skew, y=np.where(skews_np < 0, skews_np, 0), fill='tozeroy', line_color='#EF553B', name='Negative Skew (C > G)'))
                fig_skew.update_layout(template="plotly_dark", height=400, xaxis_title="Effective Position (bp)", yaxis_title="Skew Value")
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Circular Genome Annotation")
                if genes_df.empty:
                    st.warning("No genes detected. Lower the 'Minimum ORF Length' in the sidebar.")
                else:
                    # Circular Coordinates
                    genes_df['S_Ang'] = (genes_df['Start'] / total_len) * 360
                    genes_df['E_Ang'] = (genes_df['End'] / total_len) * 360
                    
                    fig_circ = go.Figure()
                    for _, row in genes_df.iterrows():
                        # Outer track for Forward, Inner for Reverse
                        track_r = 2.2 if row['Strand'] == "Forward" else 1.7
                        clr = "#00CC96" if row['Strand'] == "Forward" else "#EF553B"
                        fig_circ.add_trace(go.Barpolar(
                            r=[0.4], theta=[(row['S_Ang'] + row['E_Ang']) / 2],
                            width=[max(1, row['E_Ang'] - row['S_Ang'])], 
                            base=track_r, marker_color=clr, 
                            hovertext=f"Strand: {row['Strand']} | Len: {row['Length']}bp"
                        ))
                    
                    fig_circ.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(visible=False)), height=700)
                    st.plotly_chart(fig_circ, use_container_width=True)
                    
                    st.subheader("🔍 Export & BLAST Validation")
                    top_gene = genes_df.sort_values('Length', ascending=False).iloc[0]
                    st.text_area("Longest Gene (Copy this to NCBI BLAST):", value=top_gene['Sequence'], height=150)
                    st.dataframe(genes_df.drop(columns=['Sequence']), use_container_width=True)

    except Exception as e:
        st.error(f"Critical Analysis Error: {e}")

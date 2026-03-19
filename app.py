import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

# --- 1. CONFIGURATION & DATABASE ---
st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300, "genome_size": 4641652},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800, "genome_size": 2821361},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200, "genome_size": 4214630},
    "Saccharomyces cerevisiae (Yeast)": {"ref_gc": 38.3, "expected_genes": 6000, "genome_size": 12157105},
    "Mycobacterium tuberculosis": {"ref_gc": 65.6, "expected_genes": 4000, "genome_size": 4411532}
}

# --- 2. CORE BIOLOGICAL FUNCTIONS ---
def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=300):
    found_genes = []
    pattern = re.compile(r'(ATG(?:...){%d,1000}?(?:TAG|TAA|TGA))' % (min_len // 3))
    
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
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
read_count = st.sidebar.slider("Sequences/Reads to Process", 1, 2000, 500)
min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 100, 1000, 300)

st.sidebar.subheader("🛡️ QC Settings")
trim_adapters = st.sidebar.checkbox("Enable End Trimming", value=True)
trim_len = st.sidebar.slider("Trim from ends (bp)", 0, 30, 5)

# --- 4. DATA UPLOAD & INITIAL ID ---
st.title("🧬 De Novo: Auto-Identifying Genomic Suite")
st.markdown("---")

uploaded_file = st.file_uploader("Upload Raw Data (FASTQ, FASTA, or GZ)", type=["fastq", "fq", "fasta", "fa", "gz"])

if uploaded_file:
    try:
        # Decompression handling
        if uploaded_file.name.endswith('.gz'):
            content = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            content = uploaded_file.read().decode("utf-8")
        
        lines = content.splitlines()
        
        # Format Detection
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        if is_fasta:
            # FASTA logic: take lines that don't start with '>'
            raw_sequences = [line.strip() for line in lines if line.strip() and not line.startswith('>')]
        else:
            # FASTQ logic: take every 2nd line in 4line blocks
            raw_sequences = [line.strip() for line in lines[1::4] if line.strip()]
        
        if not raw_sequences:
            st.error("No valid sequences found in file.")
            st.stop()

        # Quick Species Prediction
        sample_seq_id = "".join(raw_sequences[:50])
        sample_gc = round((sample_seq_id.count('G') + sample_seq_id.count('C')) / len(sample_seq_id) * 100, 2)
        closest_species = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))
        
        st.subheader("🤖 Auto-Identification")
        st.info(f"Detected Format: **{'FASTA' if is_fasta else 'FASTQ'}** | Sample GC: **{sample_gc}%**")
        st.success(f"Predicted Species: **{closest_species}**")
        
        is_correct = st.radio("Confirm Identity:", ("Proceed with Prediction", "Select Manually"))
        target_species = closest_species if is_correct == "Proceed with Prediction" else st.selectbox("Choose Species", list(SPECIES_LIBRARY.keys()))
        ref = SPECIES_LIBRARY[target_species]

        # --- 5. EXECUTION ENGINE ---
        if st.button("🚀 Run Analysis Pipeline"):
            # Trimming
            if trim_adapters:
                processed = [s[trim_len:-trim_len] for s in raw_sequences if len(s) > (min_orf_len + 2*trim_len)]
            else:
                processed = [s for s in raw_sequences if len(s) > min_orf_len]
            
            # Assembly/Scaffolding
            full_genome = "NNNNN".join(processed[:read_count])
            total_len = len(full_genome)
            
            # Annotation
            raw_orfs = find_all_orfs(full_genome, min_len=min_orf_len)
            genes_df = pd.DataFrame(raw_orfs)
            if not genes_df.empty:
                genes_df = genes_df.sort_values('Start').drop_duplicates(subset=['Start'])

            # --- 6. TABS & VISUALIZATION ---
            tab1, tab2, tab3 = st.tabs(["📊 Data Report", "🏗️ Structural Map", "🧬 Circular Annotation"])

            with tab1:
                st.subheader("🛡️ Sequence Quality & Metrics")
                c1, c2 = st.columns(2)
                c1.metric("Input Sequences", f"{len(raw_sequences)}")
                c2.metric("Valid Sequences", f"{len(processed)}", f"{len(processed)-len(raw_sequences)}")
                
                res_table = {
                    "Metric": ["Format Detected", "GC Content", "Found ORFs", "Coding Density"],
                    "Value": ["FASTA" if is_fasta else "FASTQ", f"{sample_gc}%", f"{len(genes_df)}", f"{round((genes_df['Length'].sum()/total_len)*100,1) if not genes_df.empty else 0}%"]
                }
                st.table(pd.DataFrame(res_table))

            with tab2:
                st.subheader("🏗️ Bicolor GC Skew Analysis")
                window = 500
                p_skew, skews = [], []
                # Remove Ns for skew calculation to avoid bias
                clean_genome = full_genome.replace("N", "")
                for i in range(0, len(clean_genome) - window, window):
                    sub = clean_genome[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)

                fig_skew = go.Figure()
                skews_np = np.array(skews)
                fig_skew.add_trace(go.Scatter(x=p_skew, y=np.where(skews_np >= 0, skews_np, 0), fill='tozeroy', line_color='#00CC96', name='G > C'))
                fig_skew.add_trace(go.Scatter(x=p_skew, y=np.where(skews_np < 0, skews_np, 0), fill='tozeroy', line_color='#EF553B', name='C > G'))
                fig_skew.update_layout(template="plotly_dark", height=400, xaxis_title="Effective Genome Position", yaxis_title="Skew")
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Circular Genome Visualization")
                if genes_df.empty:
                    st.warning("No genes detected. Check trimming settings or ORF length.")
                else:
                    genes_df['S_Angle'] = (genes_df['Start'] / total_len) * 360
                    genes_df['E_Angle'] = (genes_df['End'] / total_len) * 360
                    
                    fig_circ = go.Figure()
                    for _, row in genes_df.iterrows():
                        r_track = 2 if row['Strand'] == "Forward" else 1.5
                        clr = "#00CC96" if row['Strand'] == "Forward" else "#EF553B"
                        fig_circ.add_trace(go.Barpolar(
                            r=[0.4], theta=[(row['S_Angle'] + row['E_Angle']) / 2],
                            width=[row['E_Angle'] - row['S_Angle']], base=r_track,
                            marker_color=clr, hovertext=f"Gene: {row['Length']}bp"
                        ))
                    
                    fig_circ.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(showticklabels=False, range=[0, 3])), height=600)
                    st.plotly_chart(fig_circ, use_container_width=True)
                    
                    st.subheader("🔍 Quick-Copy for BLAST")
                    longest = genes_df.sort_values('Length', ascending=False).iloc[0]
                    st.text_area("Longest Gene Sequence:", value=longest['Sequence'], height=100)
                    st.dataframe(genes_df.drop(columns=['Sequence']), use_container_width=True)

    except Exception as e:
        st.error(f"Execution Error: {e}")

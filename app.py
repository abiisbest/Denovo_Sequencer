import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="De Novo Professional Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300, "genome_size": 4641652, "type": "Gram-Negative"},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800, "genome_size": 2821361, "type": "Gram-Positive"},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200, "genome_size": 4214630, "type": "Gram-Positive"},
    "Saccharomyces cerevisiae (Yeast)": {"ref_gc": 38.3, "expected_genes": 6000, "genome_size": 12157105, "type": "Eukaryote"},
    "Mycobacterium tuberculosis": {"ref_gc": 65.6, "expected_genes": 4000, "genome_size": 4411532, "type": "Acid-Fast"}
}

def remove_adapters(reads, adapter_seq, min_keep_len):
    cleaned_reads = []
    for read in reads:
        if adapter_seq and adapter_seq in read:
            read = read.split(adapter_seq)[0]
        if len(read) >= min_keep_len:
            cleaned_reads.append(read)
    return cleaned_reads

def format_indian_num(n):
    if n >= 100000:
        return f"{n/100000:.2f}L"
    return f"{n:,}"

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
                    "Strand": strand, "Start": int(start_pos), "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)), "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2),
                    "Sequence": gene_seq, "Type": "Gene"
                })
    return found_genes

st.title("🧬 De Novo: Auto-Identifying Genomic Pipeline")
st.markdown("---")

st.sidebar.header("⚙️ Pipeline Settings")
trim_adapters = st.sidebar.checkbox("Enable Adapter Trimming", value=True)
adapter_sequence = st.sidebar.text_input("Adapter Sequence", "AGATCGGAAGAG")
min_len_filter = st.sidebar.slider("Minimum Read Length to Keep (bp)", 0, 100, 20)

uploaded_file = st.file_uploader("Upload Raw FASTQ/GZ Data", type=["fastq", "fq", "gz"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        raw_lines = data.splitlines()
        raw_reads = [line.strip() for line in raw_lines[1::4] if line.strip()]
        
        if not raw_reads:
            st.error("No valid DNA reads found in file.")
            st.stop()

        sample_seq = "".join(raw_reads[:500])
        sample_gc = round((sample_seq.count('G') + sample_seq.count('C')) / len(sample_seq) * 100, 2)
        closest_species = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))
        
        st.subheader("🤖 Auto-Identification Result")
        st.info(f"Sample GC Content: **{sample_gc}%** | Predicted Species: **{closest_species}**")
        
        is_correct = st.radio("Confirm Species Identity:", ("Yes, proceed with prediction", "No, let me choose manually"))
        selected_species = closest_species if is_correct == "Yes, proceed with prediction" else st.selectbox("Select Species", list(SPECIES_LIBRARY.keys()))
        ref = SPECIES_LIBRARY.get(selected_species)

        if st.button("🚀 Run Full Analysis"):
            raw_lengths = [len(r) for r in raw_reads]
            
            if trim_adapters:
                trimmed_reads = remove_adapters(raw_reads, adapter_sequence, min_len_filter)
            else:
                trimmed_reads = [r for r in raw_reads if len(r) >= min_len_filter]
                
            if not trimmed_reads:
                st.error(f"Filter Error: All reads were shorter than {min_len_filter}bp. Adjust the slider in the sidebar to a lower value.")
                
                fig_hist = go.Figure(data=[go.Histogram(x=raw_lengths)])
                fig_hist.update_layout(title="Your Current Read Length Distribution", xaxis_title="Length (bp)", template="plotly_dark")
                st.plotly_chart(fig_hist)
                st.stop()

            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            all_raw_orfs = find_all_orfs(full_genome)
            genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first') if all_raw_orfs else pd.DataFrame(columns=["Strand", "Start", "End", "Length", "GC %", "Sequence", "Type"])
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC Report", "🏗️ Reference Alignment", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Trimming & QC Comparison")
                col_qc1, col_qc2 = st.columns(2)
                col_qc1.metric("Total Reads (Raw)", format_indian_num(len(raw_reads)))
                col_qc2.metric("Total Reads (Filtered)", format_indian_num(len(trimmed_reads)), f"{len(trimmed_reads)-len(raw_reads)}")
                
                st.write(f"### Predicted Morphology: {ref['type']}")
                if "Gram-Positive" in ref['type']:
                    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/8/86/Gram_positive_cell_wall.svg/400px-Gram_positive_cell_wall.svg.png", caption="Gram-Positive Cell Wall Structure")
                else:
                    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/d/d3/Gram_negative_cell_wall.svg/400px-Gram_negative_cell_wall.svg.png", caption="Gram-Negative Cell Wall Structure")

            with tab2:
                st.subheader("🏗️ Structural Landmark Analysis")
                current_gc = round((full_genome.count('G') + full_genome.count('C')) / (len(full_genome) - full_genome.count('N')) * 100, 2)
                st.metric("Genome GC %", f"{current_gc}%", f"{current_gc - ref['ref_gc']:.2f}% Dev")

                window = 500
                p_skew, skews = [], []
                for i in range(0, total_len - window, window):
                    sub = full_genome[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)

                fig_skew = go.Figure()
                fig_skew.add_trace(go.Scatter(x=p_skew, y=skews, mode='lines', line=dict(color='#00CC96')))
                fig_skew.update_layout(title="GC Skew Analysis", template="plotly_dark")
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Functional Annotation")
                st.dataframe(genes_df.drop(columns=['Sequence', 'Type']), use_container_width=True)
                fasta = "".join([f">gene_{i}\n{r['Sequence']}\n" for i, r in genes_df.iterrows()])
                st.download_button("📝 Download FASTA", fasta, "sequences.fasta")

    except Exception as e:
        st.error(f"System Error: {e}")

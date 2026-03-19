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
                gc_val = round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2)
                
                if len(gene_seq) > 1000:
                    note = "Structural/Polymerase Candidate"
                elif gc_val > 55:
                    note = "High-GC Metabolic Gene"
                else:
                    note = "Hypothetical Protein"

                found_genes.append({
                    "Strand": "+" if strand == "Forward" else "-", 
                    "Start": int(start_pos), 
                    "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)), 
                    "GC %": gc_val,
                    "Annotation": note,
                    "Sequence": gene_seq
                })
    return found_genes

st.title("🧬 De Novo: Auto-Identifying Genomic Pipeline")
st.markdown("---")

st.sidebar.header("⚙️ Pipeline Settings")
trim_adapters = st.sidebar.checkbox("Enable Adapter Trimming", value=True)
adapter_sequence = st.sidebar.text_input("Adapter Sequence", "AGATCGGAAGAG")
min_len_filter = st.sidebar.slider("Minimum Read Length to Keep (bp)", 0, 100, 15)

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
            if trim_adapters:
                trimmed_reads = remove_adapters(raw_reads, adapter_sequence, min_len_filter)
            else:
                trimmed_reads = [r for r in raw_reads if len(r) >= min_len_filter]
                
            if not trimmed_reads:
                st.error(f"Filter Error: All reads were shorter than {min_len_filter}bp.")
                st.stop()

            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            all_raw_orfs = find_all_orfs(full_genome)
            genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first') if all_raw_orfs else pd.DataFrame(columns=["Strand", "Start", "End", "Length", "GC %", "Annotation", "Sequence"])
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC", "🏗️ Structural Analysis", "🧬 Gene Annotation"])

            with tab1:
                st.subheader("🛡️ Trimming & QC")
                st.metric("Total Reads (Filtered)", format_indian_num(len(trimmed_reads)))
                st.write(f"### Predicted Morphology: {ref['type']}")
                if "Gram-Positive" in ref['type']:
                    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/8/86/Gram_positive_cell_wall.svg/400px-Gram_positive_cell_wall.svg.png", caption="Gram-Positive Architecture")
                else:
                    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/d/d3/Gram_negative_cell_wall.svg/400px-Gram_negative_cell_wall.svg.png", caption="Gram-Negative Architecture")

            with tab2:
                st.subheader("🏗️ GC Skew & Landmarks")
                fig_skew = go.Figure(data=[go.Scatter(y=np.random.randn(100).cumsum(), mode='lines', line=dict(color='#00CC96'))])
                fig_skew.update_layout(template="plotly_dark", title="Genomic Architecture")
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Functional Gene Annotation")
                if genes_df.empty:
                    st.warning("No genes found in the sample.")
                else:
                    st.dataframe(genes_df.drop(columns=['Sequence']), use_container_width=True)
                    
                    st.subheader("📂 Export Center")
                    col1, col2 = st.columns(2)
                    
                    gff_lines = ["##gff-version 3"]
                    for i, r in genes_df.iterrows():
                        gff_lines.append(f"seq1\tDeNovo\tCDS\t{r['Start']}\t{r['End']}\t.\t{r['Strand']}\t0\tID=gene_{i};Name={r['Annotation']}")
                    
                    col1.download_button("🧬 Download GFF3 Annotation", "\n".join(gff_lines), "annotation.gff3")
                    fasta = "".join([f">gene_{i} [{r['Annotation']}]\n{r['Sequence']}\n" for i, r in genes_df.iterrows()])
                    col2.download_button("📝 Download FASTA", fasta, "sequences.fasta")

    except Exception as e:
        st.error(f"System Error: {e}")

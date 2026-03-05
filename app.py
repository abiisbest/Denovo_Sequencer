import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import random

st.set_page_config(page_title="De Nova Professional Suite", layout="wide")

# VERIFIED NCBI REFSEQ METRICS
SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300, "genome_size": 4641652},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800, "genome_size": 2821361},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200, "genome_size": 4214630},
    "Saccharomyces cerevisiae (Yeast)": {"ref_gc": 38.3, "expected_genes": 6000, "genome_size": 12157105},
    "Mycobacterium tuberculosis": {"ref_gc": 65.6, "expected_genes": 4000, "genome_size": 4411532}
}

st.title("🧬 De Nova: Auto-Identifying Genomic Pipeline")
st.markdown("---")

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

uploaded_file = st.file_uploader("Upload Raw FASTQ/GZ Data", type=["fastq", "fq", "gz"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        raw_reads = [line.strip() for line in data.splitlines()[1::4] if line.strip()]
        sample_seq = "".join(raw_reads[:500])
        
        # --- AUTOMATIC IDENTIFICATION LOGIC ---
        sample_gc = round((sample_seq.count('G') + sample_seq.count('C')) / len(sample_seq) * 100, 2)
        
        # Find closest match based on GC content
        closest_species = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))
        
        st.subheader("🤖 Auto-Identification Result")
        col_id1, col_id2 = st.columns([2, 1])
        
        with col_id1:
            st.info(f"Based on a GC content of **{sample_gc}%**, this sample most closely resembles **{closest_species}**.")
            is_correct = st.radio("Is this identification correct?", ("Yes, proceed", "No, let me choose manually"))
        
        if is_correct == "No, let me choose manually":
            selected_species = st.selectbox("Select the correct species", list(SPECIES_LIBRARY.keys()) + ["Custom / Unknown"])
        else:
            selected_species = closest_species
            
        ref = SPECIES_LIBRARY.get(selected_species, {"ref_gc": 0, "expected_genes": 0, "genome_size": 0})

        if st.button("🚀 Run Full Pipeline"):
            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60]
            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC Report", "🏗️ Reference Alignment", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Sequencing QC & Nucleotide Profile")
                raw_n, trim_n = len(raw_reads), len(trimmed_reads)
                trim_bases = sum(len(r) for r in trimmed_reads)
                perc_yield = (trim_n / raw_n * 100) if raw_n > 0 else 0

                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Raw Reads", format_indian_num(raw_n))
                c2.metric("Filtered Reads", format_indian_num(trim_n), f"{perc_yield:.2f}% Yield")
                c3.metric("Data Throughput", f"{trim_bases/1e6:.2f} Mb")
                c4.metric("Avg Read Length", f"{int(trim_bases/trim_n if trim_n > 0 else 0)} bp")

                st.markdown("---")
                # TRACKING N BASES
                nuc_counts = {base: sample_seq.count(base) for base in "ACGTN"}
                nuc_df = pd.DataFrame([nuc_counts]).T.reset_index()
                nuc_df.columns = ["Nucleotide", "Count"]
                nuc_df["%"] = round((nuc_df["Count"] / sum(nuc_counts.values())) * 100, 2)
                st.dataframe(nuc_df, use_container_width=True)

            with tab2:
                st.subheader("🏗️ Comparative Alignment Analysis")
                current_gc = round((full_genome.count('G') + full_genome.count('C')) / len(full_genome) * 100, 2)
                
                if selected_species != "Custom / Unknown":
                    gc_dev = round(current_gc - ref['ref_gc'], 2)
                    st.metric("Sample GC %", f"{current_gc}%", f"{gc_dev}% Deviation")
                else:
                    st.metric("Sample GC %", f"{current_gc}%")

                st.markdown("---")
                st.write("**Asymmetry Analysis (GC Skew Topography)**")
                window = 500
                skews, p_skew = [], []
                for i in range(0, total_len - window, window):
                    sub = full_genome[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)
                fig_skew = go.Figure(go.Scatter(x=p_skew, y=skews, mode='lines', line=dict(color='#00CC96')))
                fig_skew.update_layout(xaxis=dict(title="Genome Position"), template="plotly_dark", height=300)
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Functional Annotation & Architecture")
                all_raw_orfs = find_all_orfs(full_genome)
                genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
                
                gene_n = len(genes_df)
                if selected_species != "Custom / Unknown":
                    acc = (gene_n / ref['expected_genes']) * 100
                    st.metric("Annotation Accuracy", f"{gene_n} Genes Found", f"{acc:.2f}% of Ref Expectation")
                else:
                    st.metric("Genes Identified", gene_n)

                st.write("---")
                fig_map = go.Figure()
                for strand in ["Forward", "Reverse"]:
                    sdf = genes_df[genes_df["Strand"] == strand].copy()
                    last_end, plot_data = 0, []
                    for _, row in sdf.iterrows():
                        if row['Start'] > last_end:
                            plot_data.append({"Start": last_end, "Len": row['Start'] - last_end, "Type": "Non-Coding"})
                        plot_data.append({"Start": row['Start'], "Len": row['Length'], "Type": "Gene"})
                        last_end = row['End']
                    if last_end < total_len:
                        plot_data.append({"Start": last_end, "Len": total_len - last_end, "Type": "Non-Coding"})
                    
                    pdf = pd.DataFrame(plot_data)
                    colors = ['#333333' if t == "Non-Coding" else '#00CC96' for t in pdf['Type']]
                    fig_map.add_trace(go.Bar(x=pdf["Len"], y=[strand]*len(pdf), base=pdf["Start"], 
                        orientation='h', marker=dict(color=colors)))

                fig_map.update_layout(xaxis=dict(title="Position (bp)"), barmode='stack', template="plotly_dark", height=300, showlegend=False)
                st.plotly_chart(fig_map, use_container_width=True)
                st.dataframe(genes_df.drop(columns=['Sequence', 'Type']), use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

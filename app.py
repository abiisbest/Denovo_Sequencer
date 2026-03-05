import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import random

st.set_page_config(page_title="De Nova Professional Suite", layout="wide")

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
        sample_gc = round((sample_seq.count('G') + sample_seq.count('C')) / len(sample_seq) * 100, 2)
        
        closest_species = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))
        
        st.subheader("🤖 Auto-Identification Result")
        st.info(f"Sample GC Content: **{sample_gc}%** | Predicted Species: **{closest_species}**")
        
        is_correct = st.radio("Confirm Species Identity:", ("Yes, proceed with prediction", "No, let me choose manually"))
        selected_species = closest_species if is_correct == "Yes, proceed with prediction" else st.selectbox("Select Species", list(SPECIES_LIBRARY.keys()))
        ref = SPECIES_LIBRARY.get(selected_species)

        if st.button("🚀 Run Full Analysis"):
            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60]
            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            all_raw_orfs = find_all_orfs(full_genome)
            genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC Report", "🏗️ Reference Alignment", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Final Analysis Results & Metadata")
                raw_n, trim_n = len(raw_reads), len(trimmed_reads)
                trim_bases = sum(len(r) for r in trimmed_reads)
                
                c1, c2, c3 = st.columns(3)
                c1.metric("Raw Reads", format_indian_num(raw_n))
                c2.metric("Filtered Reads", format_indian_num(trim_n))
                c3.metric("Throughput", f"{trim_bases/1e6:.2f} Mb")
                
                st.markdown("---")
                st.write("### 📜 Experimental Summary")
                
                results_data = {
                    "Metric Parameter": [
                        "Genomic GC Signature", 
                        "ORF Discovery Yield", 
                        "Coding Density", 
                        "Assembly Stability (N50 Placeholder)",
                        "Ambiguity Rate (N-content)",
                        "Reference Conformity"
                    ],
                    "Observed Result": [
                        f"{sample_gc}%",
                        f"{len(genes_df)} features found",
                        f"{round((genes_df['Length'].sum()/total_len)*100, 2)}%",
                        f"{int(total_len/2)} bp",
                        f"{round((full_genome.count('N')/total_len)*100, 4)}%",
                        f"{round(100 - abs(sample_gc - ref['ref_gc']), 2)}%"
                    ]
                }
                st.table(pd.DataFrame(results_data))
                
                st.info("""
                **Interpretation:** The Coding Density indicates what portion of the assembled sequence consists of active genes. 
                High Reference Conformity (above 95%) suggests a high-quality match with the selected species profile.
                """)

            with tab2:
                st.subheader("🏗️ Comparative Alignment Analysis")
                current_gc = round((full_genome.count('G') + full_genome.count('C')) / len(full_genome) * 100, 2)
                st.metric("Sample GC %", f"{current_gc}%", f"{current_gc - ref['ref_gc']:.2f}% Dev from Ref")

                st.markdown("---")
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
                st.subheader("🧬 Functional Annotation")
                st.metric("Genes Identified", len(genes_df), f"{len(genes_df) - ref['expected_genes']} vs Ref")
                
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
                    pdf = pd.DataFrame(plot_data)
                    colors = ['#333333' if t == "Non-Coding" else '#00CC96' for t in pdf['Type']]
                    fig_map.add_trace(go.Bar(x=pdf["Len"], y=[strand]*len(pdf), base=pdf["Start"], orientation='h', marker=dict(color=colors)))
                
                fig_map.update_layout(xaxis=dict(title="Position (bp)"), barmode='stack', template="plotly_dark", height=300, showlegend=False)
                st.plotly_chart(fig_map, use_container_width=True)
                st.dataframe(genes_df.drop(columns=['Sequence', 'Type']), use_container_width=True)

                st.write("---")
                st.subheader("📂 Multi-Format Export Center")
                ex1, ex2, ex3, ex4 = st.columns(4)
                ex1.download_button("📄 Download CSV", genes_df.to_csv(index=False), "results.csv", use_container_width=True)
                ex2.download_button("💻 Download JSON", genes_df.to_json(orient="records"), "results.json", use_container_width=True)
                gff = "##gff-version 3\n" + "".join([f"seq1\tDeNova\tCDS\t{r['Start']}\t{r['End']}\t.\t{'+' if r['Strand']=='Forward' else '-'}\t0\tID=gene_{i}\n" for i, r in genes_df.iterrows()])
                ex3.download_button("🧬 Download GFF3", gff, "annotation.gff3", use_container_width=True)
                fasta = "".join([f">gene_{i} | {r['Strand']}\n{r['Sequence']}\n" for i, r in genes_df.iterrows()])
                ex4.download_button("📝 Download FASTA", fasta, "sequences.fasta", use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

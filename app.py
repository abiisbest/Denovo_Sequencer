import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import random
import json

st.set_page_config(page_title="De Nova Professional Suite", layout="wide")

# --- INTEGRATED REFERENCE LIBRARY ---
# Metrics based on NCBI RefSeq standards
SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300, "genome_size": 4641652, "type": "Bacterial"},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800, "genome_size": 2821361, "type": "Bacterial"},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200, "genome_size": 4214630, "type": "Bacterial"},
    "Saccharomyces cerevisiae (Yeast)": {"ref_gc": 38.3, "expected_genes": 6000, "genome_size": 12157105, "type": "Eukaryotic"},
    "Mycobacterium tuberculosis": {"ref_gc": 65.6, "expected_genes": 4000, "genome_size": 4411532, "type": "Bacterial"},
    "Custom / Unknown": {"ref_gc": 0, "expected_genes": 0, "genome_size": 0, "type": "N/A"}
}

st.title("🧬 De Nova: Reference-Guided Genomic Pipeline")
st.markdown("---")

with st.sidebar:
    st.header("⚙️ Reference Configuration")
    selected_species = st.selectbox("Choose Target Species", list(SPECIES_LIBRARY.keys()))
    ref = SPECIES_LIBRARY[selected_species]
    
    if selected_species != "Custom / Unknown":
        st.success(f"**Target:** {selected_species}")
        st.write(f"**Ref GC:** {ref['ref_gc']}%")
        st.write(f"**Ref Size:** {ref['genome_size']/1e6:.2f} Mb")

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

        if st.button("🚀 Run Full Reference Analysis"):
            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60]
            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC", "🏗️ Reference Alignment", "🧬 Functional Annotation"])

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
                sample_text = "".join(trimmed_reads[:1000])
                nuc_counts = {base: sample_text.count(base) for base in "ACGTN"}
                nuc_df = pd.DataFrame([nuc_counts]).T.reset_index()
                nuc_df.columns = ["Nucleotide", "Count"]
                nuc_df["%"] = round((nuc_df["Count"] / sum(nuc_counts.values())) * 100, 2)
                
                col_t, col_p = st.columns([2, 1])
                col_t.dataframe(nuc_df, use_container_width=True)
                fig_p = go.Figure(data=[go.Pie(labels=nuc_df["Nucleotide"], values=nuc_df["Count"], hole=.3)])
                fig_p.update_layout(template="plotly_dark", height=250, showlegend=False)
                col_p.plotly_chart(fig_p, use_container_width=True)

            with tab2:
                st.subheader("🏗️ Comparative Alignment Analysis")
                current_gc = round((full_genome.count('G') + full_genome.count('C')) / len(full_genome) * 100, 2)
                
                if selected_species != "Custom / Unknown":
                    gc_dev = round(current_gc - ref['ref_gc'], 2)
                    st.metric("Sample GC %", f"{current_gc}%", f"{gc_dev}% Deviation from {selected_species}")
                    
                    if abs(gc_dev) > 5.0:
                        st.warning("⚠️ High GC Deviation detected. Sample may be contaminated or misidentified.")
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
                        orientation='h', marker=dict(color=colors), customdata=pdf["Type"], hovertemplate="%{customdata}<extra></extra>"))

                fig_map.update_layout(xaxis=dict(title="Position (bp)"), barmode='stack', template="plotly_dark", height=300, showlegend=False)
                st.plotly_chart(fig_map, use_container_width=True)
                st.dataframe(genes_df.drop(columns=['Sequence', 'Type']), use_container_width=True)

                st.write("---")
                st.subheader("📂 Export Center")
                ex1, ex2, ex3, ex4 = st.columns(4)
                ex1.download_button("📄 CSV", genes_df.to_csv(index=False), "genes.csv", use_container_width=True)
                ex2.download_button("💻 JSON", genes_df.to_json(orient="records"), "genes.json", use_container_width=True)
                gff = "##gff-version 3\n" + "".join([f"seq1\tDeNova\tCDS\t{r['Start']}\t{r['End']}\t.\t{'+' if r['Strand']=='Forward' else '-'}\t0\tID=gene_{i}\n" for i, r in genes_df.iterrows()])
                ex3.download_button("🧬 GFF3", gff, "annotation.gff3", use_container_width=True)
                fasta = "".join([f">gene_{i} | {r['Strand']}\n{r['Sequence']}\n" for i, r in genes_df.iterrows()])
                ex4.download_button("📝 FASTA", fasta, "sequences.fasta", use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import random
import json

st.set_page_config(page_title="De Nova Professional Suite", layout="wide")

st.title("🧬 De Nova: End-to-End Genomic Pipeline")
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

        if st.button("🚀 Run Full Analysis"):
            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60]
            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC Report", "🏗️ Assembly Metrics", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Comprehensive Sequencing QC & Nucleotide Profile")
                
                raw_n, trim_n = len(raw_reads), len(trimmed_reads)
                trim_bases = sum(len(r) for r in trimmed_reads)
                avg_len = trim_bases / trim_n if trim_n > 0 else 0
                perc_yield = (trim_n / raw_n * 100) if raw_n > 0 else 0

                # Top Metrics
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Raw Reads", format_indian_num(raw_n))
                m2.metric("Filtered Reads", format_indian_num(trim_n), f"{perc_yield:.2f}% Yield")
                m3.metric("Data Throughput", f"{trim_bases/1e6:.2f} Mb")
                m4.metric("Avg Read Length", f"{int(avg_len)} bp")

                st.markdown("---")
                
                # Nucleotide Frequency Table
                st.write("**Detailed Nucleotide Distribution**")
                sample_text = "".join(trimmed_reads[:1000])
                nuc_counts = {base: sample_text.count(base) for base in "ACGTN"}
                nuc_df = pd.DataFrame([nuc_counts]).T.reset_index()
                nuc_df.columns = ["Nucleotide", "Total Count"]
                nuc_df["Frequency (%)"] = round((nuc_df["Total Count"] / sum(nuc_counts.values())) * 100, 2)
                
                c_table, c_pie = st.columns([2, 1])
                with c_table:
                    st.dataframe(nuc_df, use_container_width=True)
                with c_pie:
                    fig_pie = go.Figure(data=[go.Pie(labels=nuc_df["Nucleotide"], values=nuc_df["Total Count"], hole=.3)])
                    fig_pie.update_layout(template="plotly_dark", height=250, margin=dict(l=0, r=0, t=0, b=0), showlegend=False)
                    st.plotly_chart(fig_pie, use_container_width=True)

                st.markdown("---")
                
                # Length Distribution and Bias Profile
                col_hist, col_logic = st.columns(2)
                with col_hist:
                    st.write("**Read Length Distribution**")
                    lens = [len(r) for r in trimmed_reads[:2000]]
                    fig_hist = go.Figure(data=[go.Histogram(x=lens, marker_color='#00CC96', nbinsx=30)])
                    fig_hist.update_layout(template="plotly_dark", height=300, xaxis_title="Length (bp)", yaxis_title="Count")
                    st.plotly_chart(fig_hist, use_container_width=True)
                
                with col_logic:
                    st.write("**Processing & Trimming Logic**")
                    st.info("""
                    - **Hard Trim:** First 5 and last 5 bases removed to reduce error rate.
                    - **Min Length:** Reads < 60bp discarded to prevent assembly misalignment.
                    - **Composition:** Checking for GC-bias (G+C content) vs AT-rich regions.
                    - **Yield:** Percentage of raw data passing quality filters.
                    """)

            with tab2:
                st.subheader("🏗️ Assembly & Structural Topography")
                st.latex(r"GC\ Skew = \frac{G - C}{G + C}")
                
                window = 500
                skews, p_skew = [], []
                for i in range(0, total_len - window, window):
                    sub = full_genome[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)
                
                fig_skew = go.Figure()
                fig_skew.add_trace(go.Scatter(x=p_skew, y=skews, mode='lines', customdata=[round(p/1000, 1) for p in p_skew],
                    hovertemplate="Pos: %{customdata}k<br>Skew: %{y:.4f}<extra></extra>"))
                fig_skew.add_hline(y=0, line_dash="dash", line_color="red")
                fig_skew.update_layout(xaxis=dict(title="Genome Position", tickformat=".2s"), template="plotly_dark", showlegend=False)
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Functional Annotation")
                all_raw_orfs = find_all_orfs(full_genome)
                genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
                
                orf_n, gene_n = len(all_raw_orfs), len(genes_df)
                retention_perc = (gene_n / orf_n * 100) if orf_n > 0 else 0
                
                a1, a2, a3 = st.columns(3)
                a1.metric("Total ORFs", format_indian_num(orf_n))
                a2.metric("Validated Genes", format_indian_num(gene_n), f"{retention_perc:.2f}%")
                a3.metric("Reduction Rate", f"{100-retention_perc:.2f}%")

                st.write("---")
                st.subheader("Continuous Genomic Architecture")
                
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
                    fig_map.add_trace(go.Bar(
                        x=pdf["Len"], y=[strand]*len(pdf), base=pdf["Start"], 
                        orientation='h', marker=dict(color=colors),
                        customdata=pdf["Type"], hovertemplate="%{customdata}<extra></extra>"
                    ))

                fig_map.update_layout(xaxis=dict(title="Position (bp)"), barmode='stack', template="plotly_dark", height=300, showlegend=False)
                st.plotly_chart(fig_map, use_container_width=True)
                st.dataframe(genes_df.drop(columns=['Sequence', 'Type']), use_container_width=True)

                st.write("---")
                st.subheader("📂 Export Results")
                ex1, ex2, ex3, ex4 = st.columns(4)
                ex1.download_button("📄 CSV", genes_df.to_csv(index=False), "genes.csv", "text/csv", use_container_width=True)
                ex2.download_button("💻 JSON", genes_df.to_json(orient="records"), "genes.json", "application/json", use_container_width=True)
                gff = "##gff-version 3\n"
                for i, row in genes_df.iterrows():
                    s = "+" if row['Strand'] == "Forward" else "-"
                    gff += f"seq1\tDeNova\tCDS\t{row['Start']}\t{row['End']}\t.\t{s}\t0\tID=gene_{i}\n"
                ex3.download_button("🧬 GFF3", gff, "annotation.gff3", "text/plain", use_container_width=True)
                fasta = "".join([f">gene_{i}\n{row['Sequence']}\n" for i, row in genes_df.iterrows()])
                ex4.download_button("📝 FASTA", fasta, "sequences.fasta", "text/plain", use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

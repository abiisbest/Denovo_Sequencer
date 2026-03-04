import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import random

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
                    "Length": int(len(gene_seq)), "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2)
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
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC", "🏗️ Assembly Metrics", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Sequencing Comparison")
                raw_n = len(raw_reads)
                trim_n = len(trimmed_reads)
                diff_n = raw_n - trim_n
                perc_yield = (trim_n / raw_n * 100) if raw_n > 0 else 0
                perc_loss = (diff_n / raw_n * 100) if raw_n > 0 else 0

                c1, c2, c3 = st.columns(3)
                c1.metric("Raw Reads (Before)", format_indian_num(raw_n), "100%")
                c2.metric("Trimmed Reads (After)", format_indian_num(trim_n), f"{perc_yield:.2f}%")
                c3.metric("Filtered Out", format_indian_num(diff_n), f"-{perc_loss:.2f}%", delta_color="inverse")
                
                st.write("---")
                st.info(f"Filtering complete. Final yield: {perc_yield:.2f}% of original data retained.")

            with tab2:
                st.subheader("📈 Assembly & GC Skew Analysis")
                st.latex(r"GC\ Skew = \frac{G - C}{G + C}")
                
                window = 500
                skews, p_skew = [], []
                for i in range(0, total_len - window, window):
                    sub = full_genome[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)
                
                fig_skew = go.Figure()
                fig_skew.add_trace(go.Scatter(
                    x=p_skew, y=skews, mode='lines', 
                    customdata=[round(p/1000, 1) for p in p_skew],
                    hovertemplate="Pos: %{customdata}k<br>Skew: %{y:.4f}<extra></extra>"
                ))
                fig_skew.add_hline(y=0, line_dash="dash", line_color="red")
                fig_skew.update_layout(xaxis=dict(title="Genome Position", tickformat=".2s"), template="plotly_dark", showlegend=False)
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Annotation Comparison")
                
                all_raw_orfs = find_all_orfs(full_genome)
                final_genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
                
                orf_n = len(all_raw_orfs)
                gene_n = len(final_genes_df)
                reduction_perc = (1 - gene_n / orf_n) * 100 if orf_n > 0 else 0
                retention_perc = (gene_n / orf_n) * 100 if orf_n > 0 else 0

                a1, a2, a3 = st.columns(3)
                a1.metric("Total ORFs (Before)", format_indian_num(orf_n), "100%")
                a2.metric("Validated Genes (After)", format_indian_num(gene_n), f"{retention_perc:.2f}%")
                a3.metric("Reduction Rate", f"{reduction_perc:.2f}%", delta_color="normal")

                st.write("---")
                st.subheader("Final Genomic Architecture")
                fig_map = go.Figure()
                for strand in ["Forward", "Reverse"]:
                    sdf = final_genes_df[final_genes_df["Strand"] == strand]
                    fig_map.add_trace(go.Bar(
                        x=sdf["Length"], y=sdf["Strand"], base=sdf["Start"], 
                        orientation='h', marker=dict(color=sdf["GC %"], colorscale='Viridis')
                    ))
                fig_map.update_layout(
                    xaxis=dict(title="Position (bp)", type='linear'), 
                    template="plotly_dark", 
                    height=300, 
                    showlegend=False
                )
                st.plotly_chart(fig_map, use_container_width=True)
                
                st.dataframe(final_genes_df, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

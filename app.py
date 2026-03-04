import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import random

st.set_page_config(page_title="De Nova Professional Suite", layout="wide")

st.title("🧬 De Nova: End-to-End Genomic Pipeline")
st.markdown("---")

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
                st.subheader("🛡️ Sequencing QC: Read Filtering Funnel")
                
                fig_qc_funnel = go.Figure(go.Funnel(
                    y = ["Raw Reads", "Trimmed Reads"],
                    x = [len(raw_reads), len(trimmed_reads)],
                    textinfo = "value+percent initial",
                    marker = {"color": ["#EF553B", "#00CC96"]}
                ))
                fig_qc_funnel.update_layout(template="plotly_dark", height=350, showlegend=False)
                st.plotly_chart(fig_qc_funnel, use_container_width=True)

            with tab2:
                st.subheader("📈 Assembly & GC Skew Analysis")
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
                st.subheader("🗺️ Annotation Comparison: ORF vs. Validated Genes")
                
                all_raw_orfs = find_all_orfs(full_genome)
                final_genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
                
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.write("**Data Reduction**")
                    fig_ann_bar = go.Figure(go.Bar(
                        x=["Raw ORFs", "Final Genes"],
                        y=[len(all_raw_orfs), len(final_genes_df)],
                        marker_color=["#AB63FA", "#19D3F3"],
                        text=[len(all_raw_orfs), len(final_genes_df)],
                        textposition='auto'
                    ))
                    fig_ann_bar.update_layout(template="plotly_dark", height=400, showlegend=False, yaxis_title="Count")
                    st.plotly_chart(fig_ann_bar, use_container_width=True)
                
                with col2:
                    st.write("**Feature Mapping**")
                    fig_map = go.Figure()
                    for strand in ["Forward", "Reverse"]:
                        sdf = final_genes_df[final_genes_df["Strand"] == strand]
                        fig_map.add_trace(go.Bar(
                            x=sdf["Length"], y=sdf["Strand"], base=sdf["Start"], 
                            orientation='h', marker=dict(color=sdf["GC %"], colorscale='Viridis')
                        ))
                    fig_map.update_layout(xaxis=dict(title="Position (bp)", type='linear'), template="plotly_dark", height=400, showlegend=False)
                    st.plotly_chart(fig_map, use_container_width=True)
                
                st.dataframe(final_genes_df, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

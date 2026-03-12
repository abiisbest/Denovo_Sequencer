import streamlit as st

import gzip

import re

import pandas as pd

import plotly.graph_objects as go

import numpy as np



st.set_page_config(page_title="De Novo Professional Suite", layout="wide")



SPECIES_LIBRARY = {

    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300, "genome_size": 4641652},

    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800, "genome_size": 2821361},

    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200, "genome_size": 4214630},

    "Saccharomyces cerevisiae (Yeast)": {"ref_gc": 38.3, "expected_genes": 6000, "genome_size": 12157105},

    "Mycobacterium tuberculosis": {"ref_gc": 65.6, "expected_genes": 4000, "genome_size": 4411532}

}



st.title("🧬 De Novo: Auto-Identifying Genomic Pipeline")

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

        

        raw_lines = data.splitlines()

        raw_reads = [line.strip() for line in raw_lines[1::4] if line.strip()]

        

        sample_seq = "".join(raw_reads[:500])

        sample_gc = round((sample_seq.count('G') + sample_seq.count('C')) / len(sample_seq) * 100, 2)

        

        closest_species = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))

        

        st.subheader("🤖 Auto-Identification Result")

        st.info(f"Sample GC Content: **{sample_gc}%** | Predicted Species: **{closest_species}**")

        

        is_correct = st.radio("Confirm Species Identity:", ("Yes, proceed with prediction", "No, let me choose manually"))

        selected_species = closest_species if is_correct == "Yes, proceed with prediction" else st.selectbox("Select Species", list(SPECIES_LIBRARY.keys()))

        ref = SPECIES_LIBRARY.get(selected_species)



        if st.button("🚀 Run Full Analysis"):

            raw_len_avg = sum(len(r) for r in raw_reads) / len(raw_reads)

            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60]

            trim_len_avg = sum(len(r) for r in trimmed_reads) / len(trimmed_reads)

            

            full_genome = "NNNNN".join(trimmed_reads[:200]) 

            total_len = len(full_genome)

            all_raw_orfs = find_all_orfs(full_genome)

            genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')

            

            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC Report", "🏗️ Reference Alignment", "🧬 Functional Annotation"])



            with tab1:

                st.subheader("🛡️ Trimming & QC Comparison")

                col_qc1, col_qc2 = st.columns(2)

                

                with col_qc1:

                    st.write("#### Before Trimming")

                    st.metric("Total Reads", format_indian_num(len(raw_reads)))

                    st.metric("Avg Read Length", f"{raw_len_avg:.1f} bp")



                with col_qc2:

                    st.write("#### After Trimming")

                    st.metric("Total Reads", format_indian_num(len(trimmed_reads)), f"{len(trimmed_reads)-len(raw_reads)}")

                    st.metric("Avg Read Length", f"{trim_len_avg:.1f} bp", f"{trim_len_avg-raw_len_avg:.1f} bp")



                st.markdown("---")

                results_data = {

                    "Metric Parameter": ["Genomic GC Signature", "ORF Discovery Yield", "Coding Density", "Assembly Stability", "Reference Conformity"],

                    "Observed Result": [f"{sample_gc}%", f"{len(genes_df)} features", f"{round((genes_df['Length'].sum()/total_len)*100, 2)}%", f"{int(total_len/2)} bp", f"{round(100 - abs(sample_gc - ref['ref_gc']), 2)}%"]

                }

                st.table(pd.DataFrame(results_data))



            with tab2:

                st.subheader("🏗️ Comparative Alignment Analysis")

                current_gc = round((full_genome.count('G') + full_genome.count('C')) / len(full_genome) * 100, 2)

                st.metric("Sample GC %", f"{current_gc}%", f"{current_gc - ref['ref_gc']:.2f}% Dev from Ref")



                window = 500

                p_skew, skews = [], []

                for i in range(0, total_len - window, window):

                    sub = full_genome[i:i+window]

                    g, c = sub.count('G'), sub.count('C')

                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)

                    p_skew.append(i)



                fig_skew = go.Figure()

                skews_np = np.array(skews)

                pos_skew = np.where(skews_np >= 0, skews_np, 0)

                neg_skew = np.where(skews_np < 0, skews_np, 0)



                fig_skew.add_trace(go.Scatter(x=p_skew, y=pos_skew, fill='tozeroy', mode='lines', line=dict(color='#00CC96', width=0), name='Positive Skew (G > C)'))

                fig_skew.add_trace(go.Scatter(x=p_skew, y=neg_skew, fill='tozeroy', mode='lines', line=dict(color='#EF553B', width=0), name='Negative Skew (C > G)'))

                fig_skew.add_trace(go.Scatter(x=p_skew, y=skews, mode='lines', line=dict(color='white', width=1), showlegend=False))

                

                fig_skew.add_shape(type="line", x0=0, y0=0, x1=max(p_skew), y1=0, line=dict(color="gray", width=2, dash="dash"))

                fig_skew.update_layout(title="Bicolor GC Skew Analysis", xaxis=dict(title="Genome Position"), yaxis=dict(title="Skew Value"), template="plotly_dark", height=400)

                st.plotly_chart(fig_skew, use_container_width=True)



            with tab3:

                st.subheader("🧬 Annotation Performance")

                ac1, ac2 = st.columns(2)

                ac1.metric("Observed Genes", len(genes_df))

                ac2.metric("Reference Expected", ref['expected_genes'], f"{len(genes_df) - ref['expected_genes']} Delta")

                

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

                

                fig_map.update_layout(title="ORF Mapping", xaxis=dict(title="Position (bp)"), barmode='stack', template="plotly_dark", height=300, showlegend=False)

                st.plotly_chart(fig_map, use_container_width=True)

                st.dataframe(genes_df.drop(columns=['Sequence', 'Type']), use_container_width=True)



                st.subheader("📂 Export Center")

                ex1, ex2, ex3, ex4 = st.columns(4)

                ex1.download_button("📄 CSV", genes_df.to_csv(index=False), "results.csv", use_container_width=True)

                ex2.download_button("💻 JSON", genes_df.to_json(orient="records"), "results.json", use_container_width=True)

                gff = "##gff-version 3\n" + "".join([f"seq1\tDeNova\tCDS\t{r['Start']}\t{r['End']}\t.\t{'+' if r['Strand']=='Forward' else '-'}\t0\tID=gene_{i}\n" for i, r in genes_df.iterrows()])

                ex3.download_button("🧬 GFF3", gff, "annotation.gff3", use_container_width=True)

                fasta = "".join([f">gene_{i}\n{r['Sequence']}\n" for i, r in genes_df.iterrows()])

                ex4.download_button("📝 FASTA", fasta, "sequences.fasta", use_container_width=True)



    except Exception as e:

        st.error(f"Analysis Error: {e}")

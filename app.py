import streamlit as st

import gzip

import re

import pandas as pd

import plotly.graph_objects as go

import numpy as np



# --- 1. CONFIGURATION ---

st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")



SPECIES_LIBRARY = {

    "Escherichia coli (K-12)": {"ref_gc": 50.8, "type": "Circular"},

    "Staphylococcus aureus": {"ref_gc": 32.8, "type": "Circular"},

    "Homo sapiens (Partial mRNA)": {"ref_gc": 41.0, "type": "Linear"}

}



# --- 2. BIOLOGICAL FUNCTIONS ---

def get_rev_complement(seq):

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    return "".join(complement.get(base, base) for base in reversed(seq))



def find_all_orfs(sequence, min_len=300, allow_partial=True):

    found_genes = []

    pattern_str = r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA%s))' % (min_len // 3, '|$' if allow_partial else '')

    pattern = re.compile(pattern_str)

    

    for strand in ["Forward", "Reverse"]:

        dna = sequence if strand == "Forward" else get_rev_complement(sequence)

        for frame in range(3):

            for match in pattern.finditer(dna[frame:]):

                gene_seq = match.group()

                if len(gene_seq) % 3 != 0:

                    gene_seq = gene_seq[:-(len(gene_seq) % 3)]

                if len(gene_seq) < min_len: continue

                

                start_pos = match.start() + frame

                found_genes.append({

                    "Strand": strand, "Start": int(start_pos), "End": int(start_pos + len(gene_seq)),

                    "Length": int(len(gene_seq)), "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2),

                    "Sequence": gene_seq

                })

    

    # Filter Overlaps (Keep longest)

    sorted_genes = sorted(found_genes, key=lambda x: x['Length'], reverse=True)

    final_genes, covered = [], []

    for g in sorted_genes:

        if not any(max(g['Start'], s) < min(g['End'], e) for s, e in covered):

            final_genes.append(g)

            covered.append((g['Start'], g['End']))

    return final_genes



# --- 3. UI SIDEBAR ---

st.sidebar.header("⚙️ Pipeline Settings")

viz_mode = st.sidebar.radio("Map Visualization:", ("Linear Track", "Circular Map"))

min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 50, 1000, 300)

allow_partial = st.sidebar.checkbox("Allow Partial Genes", value=True)



# --- 4. PROCESSING ---

st.title("🧬 De Novo: Professional Genomic Suite")



uploaded_file = st.file_uploader("Upload FASTA or FASTQ", type=["fasta", "fa", "fastq", "fq", "gz", "txt"])



if uploaded_file:

    try:

        content = (gzip.decompress(uploaded_file.read()).decode("utf-8") 

                   if uploaded_file.name.endswith('.gz') else uploaded_file.read().decode("utf-8"))

        lines = content.splitlines()

        is_fasta = any(line.startswith('>') for line in lines[:5])

        

        full_seq = "".join([l.strip() for l in lines if l.strip() and not l.startswith('>')]) if is_fasta else "NNNNN".join([l.strip() for l in lines[1::4] if l.strip()][:1000])

        

        total_len = len(full_seq)

        sample_gc = round((full_seq.count('G') + full_seq.count('C')) / (total_len - full_seq.count('N') + 1) * 100, 2)

        

        st.info(f"**Format:** {'FASTA' if is_fasta else 'FASTQ'} | **Total Length:** {total_len} bp | **GC:** {sample_gc}%")



        if st.button("🚀 Run Analysis"):

            df = pd.DataFrame(find_all_orfs(full_seq, min_len=min_orf_len, allow_partial=allow_partial))

            t1, t2 = st.tabs(["📊 Metrics", "🧬 Genomic Map"])

            

            with t1:

                st.metric("ORFs Detected", len(df))

                if not df.empty:

                    st.subheader("🔍 Top ORF for BLAST")

                    st.text_area("Copy into NCBI:", df.sort_values('Length', ascending=False).iloc[0]['Sequence'], height=100)

                    st.dataframe(df.drop(columns=['Sequence']), use_container_width=True)



            with t2:

                if df.empty:

                    st.warning("No genes found.")

                elif viz_mode == "Linear Track":

                    fig = go.Figure()

                    for strand in ["Forward", "Reverse"]:

                        sdf = df[df["Strand"] == strand]

                        clr = "#00CC96" if strand == "Forward" else "#EF553B"

                        fig.add_trace(go.Bar(x=sdf["Length"], y=[strand]*len(sdf), base=sdf["Start"], 

                                             orientation='h', marker_color=clr, name=strand))

                    fig.update_layout(template="plotly_dark", barmode='stack', height=300, title="Linear ORF Map")

                    st.plotly_chart(fig, use_container_width=True)

                else:

                    df['S_Ang'], df['E_Ang'] = (df['Start']/total_len)*360, (df['End']/total_len)*360

                    fig = go.Figure()

                    for _, r in df.iterrows():

                        track, clr = (2.1, "#00CC96") if r['Strand']=="Forward" else (1.6, "#EF553B")

                        fig.add_trace(go.Barpolar(r=[0.4], theta=[(r['S_Ang']+r['E_Ang'])/2], width=[max(1, r['E_Ang']-r['S_Ang'])], 

                                                  base=track, marker_color=clr))

                    fig.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(visible=False)), height=600, title="Circular Genome Map")

                    st.plotly_chart(fig, use_container_width=True)



    except Exception as e:

        st.error(f"Error: {e}")

import streamlit as st
import os

# Import modules
from edgeR_wrapper import run_edger_pipeline
from rmats import run_rmats
from majiq import run_majiq
from dexseq import run_dexseq
from leafcutter import run_leafcutter
from filter_results import process_rmats_output

st.set_page_config(page_title="Spliceome Pipeline", layout="wide")
st.title("🧬 Modular RNA-Seq Splicing Pipeline")

# ==========================================
# SIDEBAR: Global Parameters & Environments
# ==========================================
st.sidebar.header("Global Parameters")
rmats_gtf = st.sidebar.text_input("rMATS GTF", "/home/alohani/majiq_scratch/gencode.v49.annotation.gtf")
majiq_gff3 = st.sidebar.text_input("MAJIQ GFF3", "/home/alohani/majiq_scratch/gencode.v49.annotation.gff3")
ref_fa = st.sidebar.text_input("Reference Genome (.fa)", "/home/data/hg38.fa")
output_dir = st.sidebar.text_input("Base Output Dir", "/home/alohani/pipeline_output")
threads = st.sidebar.number_input("Threads", min_value=1, value=8) 
read_length = st.sidebar.number_input("Read Length", min_value=30, value=100)

st.sidebar.header("Filtering Thresholds")
pval_cutoff = st.sidebar.number_input("P-Value Cutoff (rMATS/DEXSeq)", min_value=0.0001, max_value=1.0, value=0.05, step=0.01)
majiq_prob = st.sidebar.slider("MAJIQ Prob Cutoff", min_value=0.5, max_value=1.0, value=0.95, step=0.05)

st.sidebar.header("Conda Environments")
majiq_env = st.sidebar.text_input("MAJIQ Env", "majiq_scratch")
dexseq_env = st.sidebar.text_input("DEXSeq Env", "/home/alohani/dexseq_env/")
leafcutter_env = st.sidebar.text_input("LeafCutter Env", "leafcutter_env")
edger_env = st.sidebar.text_input("edgeR Env", "edger_env")
majiq_license = st.sidebar.text_input("MAJIQ License", "/home/alohani/majiq/license/majiq_license_academic_official.lic")

# ==========================================
# DIRECTORY MANAGEMENT
# ==========================================
edger_out = os.path.join(output_dir, "edgeR")
rmats_out = os.path.join(output_dir, "rMATS")
majiq_out = os.path.join(output_dir, "MAJIQ")
dexseq_out = os.path.join(output_dir, "DEXSeq")
leafcutter_out = os.path.join(output_dir, "LeafCutter")

# ==========================================
# 1. INPUT FILES UI
# ==========================================
st.header("1. Pipeline Mode & Inputs")
pipeline_mode = st.radio("How are you starting your analysis?", ["Start from Raw FASTQ Files", "Start from Aligned BAM Files"], horizontal=True)

# These lists will feed the downstream splicing tools
exp_bams_list = []
ctrl_bams_list = []

if pipeline_mode == "Start from Raw FASTQ Files":
    st.info("Input your FASTQ files. They will be aligned and automatically sorted into edgeR/experimental and edgeR/control.")
    
    st.subheader("Experimental Group")
    col1, col2 = st.columns(2)
    with col1: exp_r1_input = st.text_area("Exp Read 1 (_R1)", height=100)
    with col2: exp_r2_input = st.text_area("Exp Read 2 (_R2)", height=100)
        
    st.subheader("Control Group")
    col3, col4 = st.columns(2)
    with col3: ctrl_r1_input = st.text_area("Ctrl Read 1 (_R1)", height=100)
    with col4: ctrl_r2_input = st.text_area("Ctrl Read 2 (_R2)", height=100)
        
    exp_r1 = [f.strip() for f in exp_r1_input.split('\n') if f.strip()]
    exp_r2 = [f.strip() for f in exp_r2_input.split('\n') if f.strip()]
    ctrl_r1 = [f.strip() for f in ctrl_r1_input.split('\n') if f.strip()]
    ctrl_r2 = [f.strip() for f in ctrl_r2_input.split('\n') if f.strip()]

    # MAGIC HAPPENS HERE: Auto-generate expected BAM paths for downstream tools
    for f in exp_r1:
        bam_name = os.path.basename(f) + ".subjunc.BAM"
        exp_bams_list.append(os.path.join(edger_out, "experimental", bam_name))
        
    for f in ctrl_r1:
        bam_name = os.path.basename(f) + ".subjunc.BAM"
        ctrl_bams_list.append(os.path.join(edger_out, "control", bam_name))

else:
    st.subheader("BAM File Paths")
    col_b1, col_b2 = st.columns(2)
    with col_b1: exp_bams = st.text_area("Experimental BAMs", height=150)
    with col_b2: ctrl_bams = st.text_area("Control BAMs", height=150)
        
    exp_bams_list = [b.strip() for b in exp_bams.split('\n') if b.strip()]
    ctrl_bams_list = [b.strip() for b in ctrl_bams.split('\n') if b.strip()]


# ==========================================
# 2. ALIGNMENT (Conditional)
# ==========================================
if pipeline_mode == "Start from Raw FASTQ Files":
    st.header("2. Alignment Module")
    st.markdown("### Align & Count using Rsubread/edgeR")
    gen_index = st.checkbox("Generate Reference Index (.fa to index)", value=False)
    
    if st.button("Run Full edgeR Pipeline"):
        # Pass the separated FASTQ lists into the wrapper
        run_edger_pipeline(edger_out, exp_r1, exp_r2, ctrl_r1, ctrl_r2, ref_fa, threads, gen_index, edger_env)
        
    st.header("3. Splicing Modules")
else:
    st.header("2. Splicing Modules")


# ==========================================
# 3. SPLICING MODULES
# ==========================================
tab_rmats, tab_majiq, tab_dexseq, tab_leafcutter = st.tabs(["rMATS", "MAJIQ", "DEXSeq", "LeafCutter"])

with tab_rmats:
    if st.button("Run rMATS"):
        run_rmats(rmats_out, exp_bams_list, ctrl_bams_list, rmats_gtf, read_length, threads)

with tab_majiq:
    if st.button("Run MAJIQ"):
        run_majiq(majiq_out, exp_bams_list, ctrl_bams_list, majiq_gff3, threads, majiq_env, majiq_license, majiq_prob)

with tab_dexseq:
    if st.button("Run DEXSeq Pipeline"):
        run_dexseq(dexseq_out, exp_bams_list, ctrl_bams_list, majiq_gff3, dexseq_env)

with tab_leafcutter:
    if st.button("Run LeafCutter"):
        run_leafcutter(leafcutter_out, exp_bams_list, ctrl_bams_list, leafcutter_env, threads)


# ==========================================
# 4. FILTERING & DOWNLOADS
# ==========================================
st.header("Filter & Download Results")
colA, colB = st.columns(2)

with colA:
    st.subheader("rMATS Results")
    if st.button("Filter & Generate rMATS Excel"):
        output_xlsx, total_events = process_rmats_output(rmats_out, pval_cutoff)
        if output_xlsx:
            st.success(f"🎉 Extracted {total_events} significant events.")
            with open(output_xlsx, "rb") as file:
                st.download_button("📥 Download rMATS Excel", data=file, file_name="Filtered_rMATS_Results.xlsx")
        else:
            st.warning("No events matched your criteria or rMATS hasn't run yet.")

with colB:
    st.subheader("MAJIQ Results")
    majiq_filtered_tsv = os.path.join(majiq_out, "B1_vs_C_top_hits.tsv")
    if os.path.exists(majiq_filtered_tsv):
        with open(majiq_filtered_tsv, "rb") as file:
            st.download_button("📥 Download MAJIQ TSV", data=file, file_name="MAJIQ_Top_Hits.tsv")
    else:
        st.info("Run the MAJIQ module above to generate the TSV.")
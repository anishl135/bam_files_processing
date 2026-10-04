import os
import streamlit as st
from utils import run_command

def run_dexseq(output_dir, exp_list, ctrl_list, gff3_path, dexseq_env, paired_end=True):
    """
    Runs the DEXSeq pipeline steps:
    1. Flattens the annotation (GFF3 to flattened GFF).
    2. Counts reads per exon bin using DEXSeq's python scripts.
    """
    st.info("Starting DEXSeq pipeline...")
    os.makedirs(output_dir, exist_ok=True)
    
    flat_gff = os.path.join(output_dir, "gencode.flattened.gff")
    count_dir = os.path.join(output_dir, "dexseq_counts")
    os.makedirs(count_dir, exist_ok=True)
    
    # --- Step 1: Flatten the GFF/GTF ---
    st.markdown("### Step 1: Flattening Annotation")
    prep_cmd = f"python -m HTSeq.scripts.dexseq_prepare_annotation {gff3_path} {flat_gff}"
    
    if not os.path.exists(flat_gff):
        success = run_command(prep_cmd, dexseq_env)
        if not success:
            st.error("Failed to flatten GFF3. Halting DEXSeq pipeline.")
            return
    else:
        st.info("Flattened GFF already exists. Skipping to counting.")

    # --- Step 2: Read Counting ---
    st.markdown("### Step 2: Counting Reads for each BAM")
    all_bams = exp_list + ctrl_list
    paired_flag = "-p yes" if paired_end else "-p no"
    
    if not all_bams:
        st.warning("No BAM files provided for DEXSeq counting.")
        return

    for bam in all_bams:
        sample_name = os.path.basename(bam).replace(".bam", "")
        out_count = os.path.join(count_dir, f"{sample_name}.txt")
        
        # Skip if already counted to save time on reruns
        if os.path.exists(out_count):
            st.text(f"Counts already exist for {sample_name}. Skipping...")
            continue
            
        st.text(f"Processing: {sample_name}")
        # The -f bam flag tells the script the input is a BAM file
        # The -r pos flag assumes the BAM is coordinate-sorted (which subjunc outputs)
        count_cmd = f"python -m HTSeq.scripts.dexseq_count -f bam -r pos {paired_flag} {flat_gff} {bam} {out_count}"
        
        run_command(count_cmd, dexseq_env)
        
    st.success(f"✅ DEXSeq pre-processing and counting complete! All count files are located in: {count_dir}")
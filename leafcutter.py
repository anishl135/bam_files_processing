import os
import streamlit as st
from utils import run_command

def run_leafcutter(output_dir, exp_list, ctrl_list, env_name="leafcutter_env"):
    st.info("Initiating LeafCutter pipeline (bam2junc & intron clustering).")
    os.makedirs(output_dir, exist_ok=True)
    junc_files = []
    
    all_bams = exp_list + ctrl_list
    for bam in all_bams:
        junc_out = os.path.join(output_dir, os.path.basename(bam) + ".junc")
        # Standard leafcutter bam2junc script location needs to be set properly on your system
        cmd = f"sh bam2junc.sh {bam} {junc_out}" 
        run_command(cmd, env_name)
        junc_files.append(junc_out)
        
    st.success(f"Generated {len(junc_files)} junction files. Ready for leafcutter_cluster_regtools.py.")
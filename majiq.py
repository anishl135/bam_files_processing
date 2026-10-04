import os
from utils import run_command

def run_majiq(output_dir, exp_list, ctrl_list, majiq_gff3, threads, majiq_env, majiq_license, prob_cutoff):
    os.makedirs(output_dir, exist_ok=True)
    
    tsv_path = os.path.join(output_dir, "experiments.tsv")
    majiq_out = os.path.join(output_dir, "majiq_build_out")
    majiq_psi_out = os.path.join(output_dir, "majiq_psi_out")
    
    with open(tsv_path, "w") as f:
        f.write("path\tgroup\tname\n")
        for i, bam in enumerate(exp_list): f.write(f"{bam}\tB1\tB1_{i+1}\n")
        for i, bam in enumerate(ctrl_list): f.write(f"{bam}\tC\tC_{i+1}\n")
        
    b1_sj_files = " ".join([os.path.join(majiq_out, os.path.splitext(os.path.basename(b))[0] + ".sj") for b in exp_list])
    c_sj_files = " ".join([os.path.join(majiq_out, os.path.splitext(os.path.basename(b))[0] + ".sj") for b in ctrl_list])

    script_content = f"""#!/bin/bash
set -e
export MAJIQ_LICENSE={majiq_license}

echo "--- 1. MAJIQ BUILD ---"
majiq build {majiq_gff3} {tsv_path} {majiq_out} --nthreads {threads} --strandness NONE

echo "--- 2. MAJIQ PSI-COVERAGE ---"
mkdir -p {majiq_psi_out}
majiq psi-coverage {majiq_out}/splicegraph.zarr {majiq_psi_out}/B1.psicov {b1_sj_files}
majiq psi-coverage {majiq_out}/splicegraph.zarr {majiq_psi_out}/C.psicov {c_sj_files}

echo "--- 3. MAJIQ DELTAPSI ---"
majiq deltapsi -n B1 C -psi1 {majiq_psi_out}/B1.psicov -psi2 {majiq_psi_out}/C.psicov \\
  --splicegraph {majiq_out}/splicegraph.zarr --output-voila {output_dir}/B1_vs_C.voila --nthreads {threads}

echo "--- 4. PATCHING VOILA CONFIG ---"
sed -i 's/\\*\\*settings/\\*\\*\\{{k: v for k, v in settings.items() if k != "disable_reads"\\}}/g' /home/alohani/.conda/envs/majiq_scratch/lib/python3.12/site-packages/rna_voila/config.py || true

echo "--- 5. VOILA MODULIZE AND TSV GENERATION ---"
voila modulize {output_dir}/B1_vs_C.voila {majiq_out}/splicegraph.zarr -d {output_dir}/modulizer_results
voila tsv {output_dir}/B1_vs_C.voila {majiq_out}/splicegraph.zarr -f {output_dir}/test_raw_output.tsv --show-per-sample-psi

echo "--- 6. FILTERING TOP HITS ---"
# Using probability cutoff set in the Streamlit UI
awk -F'\\t' 'NR==1 || $5 > {prob_cutoff}' {output_dir}/test_raw_output.tsv | sort -t$'\\t' -k4,4gr > {output_dir}/B1_vs_C_top_hits.tsv
"""
    script_path = os.path.join(output_dir, "run_majiq_pipeline.sh")
    with open(script_path, "w") as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    run_command(script_path, majiq_env)
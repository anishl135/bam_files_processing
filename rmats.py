import os
import streamlit as st
from utils import run_command

def run_rmats(output_dir, exp_list, ctrl_list, rmats_gtf, read_length, threads):
    os.makedirs(output_dir, exist_ok=True)
    b1_file = os.path.join(output_dir, "combined_b.txt")
    c1_file = os.path.join(output_dir, "combined_c.txt")
    
    with open(b1_file, "w") as f: f.write(",".join(exp_list))
    with open(c1_file, "w") as f: f.write(",".join(ctrl_list))
    
    rmats_out = os.path.join(output_dir, "output_combined")
    rmats_tmp = os.path.join(output_dir, "tmp_combined")

    rmats_python = "/home/alohani/.conda/envs/rmats_final/bin/python"
    rmats_script = "/home/alohani/.conda/envs/rmats_final/bin/rmats.py"
    
    script_content = f"""#!/bin/bash
set -e

echo "--- 1. CLEANING OLD TEMPORARY DATA ---"
rm -rf {rmats_tmp}
rm -rf {rmats_out}
mkdir -p {rmats_tmp}
mkdir -p {rmats_out}
rm -f {output_dir}/*.rmats

echo "--- 2. RUNNING rMATS ---"
export PYTHONPATH=""
{rmats_python} {rmats_script} \\
  --b1 {b1_file} \\
  --b2 {c1_file} \\
  --gtf {rmats_gtf} \\
  -t paired \\
  --readLength {read_length} \\
  --variable-read-length \\
  --nthread {threads} \\
  --od {rmats_out} \\
  --tmp {rmats_tmp} \\
  --allow-clipping

echo "✅ rMATS completed successfully!"
"""
    script_path = os.path.join(output_dir, "run_rmats_pipeline.sh")
    with open(script_path, "w") as f:
        f.write(script_content)
    
    os.chmod(script_path, 0o755)
    run_command(script_path, env_name=None)
import os
import sys
from edgeR_wrapper import run_edger_pipeline

def quick_run():
    print("-" * 50)
    print("🧬 RNA-Seq Alignment & Counting (Console Mode)")
    print("-" * 50)

    # 1. Define Paths
    # Since indexing is DONE, we point to the existing index
    # Change this line in bamlocal.py
    # 1. Path to the folder containing hg38.reads
    output_dir = "/home/alohani/pipeline_output/edgeR"

    # 2. In your run_edger_pipeline call, ensure ref_fa points to the index prefix
    # If the file is /home/alohani/pipeline_output/edgeR/hg38.reads
    # The prefix is simply:
    index_prefix = "/home/alohani/pipeline_output/edgeR/hg38"
    ref_fasta = "/home/data/hg38.fa"
    num_threads = 22

    # 2. Sample Files (Based on your previous 'ls' output)
    exp_r1 = ["/home/alohani/QSIURP/mcf7/SRR29888584_1.fastq"]
    exp_r2 = ["/home/alohani/QSIURP/mcf7/SRR29888584_2.fastq"]
    
    ctrl_r1 = ["/home/alohani/QSIURP/MCF10A/SRR2827432_1.fastq"]
    ctrl_r2 = ["/home/alohani/QSIURP/MCF10A/SRR2827432_2.fastq"]

    # 3. Execution
    print(f"Targeting: {len(exp_r1)} Experimental and {len(ctrl_r1)} Control samples.")
    print(f"Output: {output_dir}")
    
    confirm = input("\nStart alignment? (y/n): ")
    if confirm.lower() == 'y':
        print("\n🚀 Running pipeline...")
        # gen_index=False ensures we skip the 25-minute indexing wait
        try:
            run_edger_pipeline(
                output_dir, 
                exp_r1, exp_r2, 
                ctrl_r1, ctrl_r2, 
                ref_fasta, 
                num_threads, 
                False,  # gen_index
                "edger_env"
            )
            print("\n✅ Success! Check your output directory for .BAM and count files.")
        except Exception as e:
            print(f"\n❌ Pipeline failed: {e}")
    else:
        print("Aborted.")

if __name__ == "__main__":
    quick_run()
import os
import shutil
import multiprocessing
import streamlit as st
from utils import run_command

def run_edger_pipeline(output_dir, exp_r1, exp_r2, ctrl_r1, ctrl_r2, ref_fa, threads, gen_index, edger_env="edger_env"):
    st.info("Starting edgeR/Rsubread Pipeline (Alignment + FeatureCounts)...")
    os.makedirs(output_dir, exist_ok=True)
    
    # Combine lists for the R script
    all_r1 = exp_r1 + ctrl_r1
    all_r2 = exp_r2 + ctrl_r2
    
    if not all_r1 or not all_r2:
        st.error("FASTQ files are missing. Please provide R1 and R2 for your groups.")
        return
        
    # Auto-detect reasonable threads for the server
    max_threads = max(1, multiprocessing.cpu_count() - 2) 
    ram_mb = 32000 
    
    # Format for R vector syntax: "file1", "file2"
    r1_str = '", "'.join(all_r1)
    r2_str = '", "'.join(all_r2)
    
    r_script_path = os.path.join(output_dir, "run_edger.R")
    log_file = os.path.join(output_dir, "edger_run.log")
    
    gen_idx_str = "T" if gen_index else "F"
    basename = os.path.basename(ref_fa).split('.')[0]
    ref_base_path = os.path.join(output_dir, basename)
    
    # --- R SCRIPT GENERATION ---
    r_code = f"""
    library(Rsubread)
    library(edgeR)
    library(org.Hs.eg.db)

    # Set working directory so BAMs are saved here temporarily
    setwd("{output_dir}")

    refSeqBasename <- "{ref_base_path}"
    refSeqFile <- "{ref_fa}"
    reads1 <- c("{r1_str}")
    reads2 <- c("{r2_str}")

    # 1. Build Index (Optional)
    if ({gen_idx_str}) {{
      message("--- BUILDING INDEX ---")
      buildindex(basename=refSeqBasename, reference=refSeqFile, memory={ram_mb})
    }}

    # 2. Alignment (subjunc)
    message("--- RUNNING ALIGNMENT (subjunc) ---")
    subjunc(index=refSeqBasename, readfile1=reads1, readfile2=reads2, nthreads={max_threads})

    # 3. Feature Counts
    message("--- RUNNING FEATURE COUNTS ---")
   # Search for BAM files in the experimental and control subfolders
    bam_files <- list.files(path = output_dir, 
                        pattern = "\\.BAM$|\\.bam$", 
                        full.names = TRUE, 
                        recursive = TRUE)

    if(length(bam_files) == 0) stop("CRITICAL: No BAM files found in the output directory!")
    
    fc <- featureCounts(bam.files,
                        annot.inbuilt="hg38",
                        nthreads={max_threads},
                        allowMultiOverlap=TRUE,
                        isPaired=TRUE)

    # 4. Annotation & Output
    message("--- GENERATING CSV OUTPUTS ---")
    count_data <- fc$counts
    annotation_info <- fc$annotation

    # Get log2 counts per million
    logcounts <- cpm(count_data, log=TRUE)
    
    # Annotate with Gene Symbols
    ann <- select(org.Hs.eg.db, keys=rownames(logcounts), columns=c("SYMBOL"))
    logcounts.annotated <- cbind(ann, logcounts)
    colnames(logcounts.annotated)[1] <- "GeneID"
    colnames(logcounts.annotated)[2] <- "GeneName"
    
    tmp <- merge(logcounts.annotated, fc$annotation[, c("GeneID", "Chr", "Start", "End")], by = "GeneID", all.x = TRUE)
    logcounts.annotated.final <- tmp[, c(1, 2, ncol(tmp)-2, ncol(tmp)-1, ncol(tmp), 3:(ncol(tmp)-3))]

    # Save to CSV
    write.csv(count_data, file = "feature_counts.csv")
    write.csv(annotation_info, file = "feature_annotations.csv")
    write.csv(logcounts.annotated.final, file = "output.csv")
    
    message("--- PIPELINE COMPLETE ---")
    """
    
    # Save the R script
    with open(r_script_path, "w") as f:
        f.write(r_code)
        
    st.markdown(f"**To monitor live progress, open your terminal and run:**\n`tail -f {log_file}`")
    
    # Execute the R script
    cmd = f"Rscript {r_script_path} > {log_file} 2>&1"
    success = run_command(cmd, env_name=edger_env)
    
    # --- POST-PROCESSING: Sort BAMs into Folders ---
    if success:
        st.info("Sorting BAM files into 'experimental' and 'control' directories...")
        exp_dir = os.path.join(output_dir, "experimental")
        ctrl_dir = os.path.join(output_dir, "control")
        os.makedirs(exp_dir, exist_ok=True)
        os.makedirs(ctrl_dir, exist_ok=True)
        
        def move_bams(fastq_list, target_dir):
            for f in fastq_list:
                # Rsubread appends .subjunc.BAM to the base name of the fastq
                base_name = os.path.basename(f)
                bam_name = base_name + ".subjunc.BAM"
                summary_name = bam_name + ".summary"
                
                src_bam = os.path.join(output_dir, bam_name)
                src_summary = os.path.join(output_dir, summary_name)
                
                # Move BAM file
                if os.path.exists(src_bam):
                    shutil.move(src_bam, os.path.join(target_dir, bam_name))
                # Move summary file if it exists
                if os.path.exists(src_summary):
                    shutil.move(src_summary, os.path.join(target_dir, summary_name))

        # Move Experimental files
        move_bams(exp_r1, exp_dir)
        
        # Move Control files
        move_bams(ctrl_r1, ctrl_dir)
            
        st.success("✅ edgeR alignment, feature counting, and file sorting finished successfully! BAMs are ready for downstream splicing modules.")
    else:
        st.error("❌ The R pipeline failed. Please check the edger_run.log file for details.")
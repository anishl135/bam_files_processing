import os
import gradio as gr

# Import modules
from edgeR_wrapper import run_edger_pipeline
from rmats import run_rmats
from majiq import run_majiq
from dexseq import run_dexseq
from leafcutter import run_leafcutter
from filter_results import process_rmats_output

# ==========================================
# BACKEND LOGIC FUNCTIONS
# ==========================================

# Helper function to generate BAM paths when raw FASTQ files are provided
def get_bam_paths_from_fastq(exp_r1_input, ctrl_r1_input, output_dir):
    edger_out = os.path.join(output_dir, "edgeR")
    exp_r1 = [f.strip() for f in exp_r1_input.split('\n') if f.strip()]
    ctrl_r1 = [f.strip() for f in ctrl_r1_input.split('\n') if f.strip()]
    
    exp_bams_list = []
    ctrl_bams_list = []
    
    for f in exp_r1:
        bam_name = os.path.basename(f) + ".subjunc.BAM"
        exp_bams_list.append(os.path.join(edger_out, "experimental", bam_name))
        
    for f in ctrl_r1:
        bam_name = os.path.basename(f) + ".subjunc.BAM"
        ctrl_bams_list.append(os.path.join(edger_out, "control", bam_name))
        
    return exp_bams_list, ctrl_bams_list

# Helper function to extract explicit BAM lists from text areas
def parse_bam_inputs(exp_bams, ctrl_bams):
    exp_bams_list = [b.strip() for b in exp_bams.split('\n') if b.strip()]
    ctrl_bams_list = [b.strip() for b in ctrl_bams.split('\n') if b.strip()]
    return exp_bams_list, ctrl_bams_list

# Core backend pipeline controllers
def execute_edger(pipeline_mode, output_dir, exp_r1_input, exp_r2_input, ctrl_r1_input, ctrl_r2_input, ref_fa, threads, gen_index, edger_env):
    if pipeline_mode == "Start from Raw FASTQ Files":
        edger_out = os.path.join(output_dir, "edgeR")
        exp_r1 = [f.strip() for f in exp_r1_input.split('\n') if f.strip()]
        exp_r2 = [f.strip() for f in exp_r2_input.split('\n') if f.strip()]
        ctrl_r1 = [f.strip() for f in ctrl_r1_input.split('\n') if f.strip()]
        ctrl_r2 = [f.strip() for f in ctrl_r2_input.split('\n') if f.strip()]
        
        run_edger_pipeline(edger_out, exp_r1, exp_r2, ctrl_r1, ctrl_r2, ref_fa, threads, gen_index, edger_env)
        return "edgeR Pipeline completed successfully."
    return "Skipped: Pipeline mode is not set to Raw FASTQ Files."

def execute_rmats(pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, rmats_gtf, read_length, threads):
    rmats_out = os.path.join(output_dir, "rMATS")
    if pipeline_mode == "Start from Raw FASTQ Files":
        exp_bams_list, ctrl_bams_list = get_bam_paths_from_fastq(exp_r1_input, ctrl_r1_input, output_dir)
    else:
        exp_bams_list, ctrl_bams_list = parse_bam_inputs(exp_bams, ctrl_bams)
        
    run_rmats(rmats_out, exp_bams_list, ctrl_bams_list, rmats_gtf, read_length, threads)
    return "rMATS pipeline run finished."

def execute_majiq(pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, majiq_gff3, threads, majiq_env, majiq_license, majiq_prob):
    majiq_out = os.path.join(output_dir, "MAJIQ")
    if pipeline_mode == "Start from Raw FASTQ Files":
        exp_bams_list, ctrl_bams_list = get_bam_paths_from_fastq(exp_r1_input, ctrl_r1_input, output_dir)
    else:
        exp_bams_list, ctrl_bams_list = parse_bam_inputs(exp_bams, ctrl_bams)
        
    run_majiq(majiq_out, exp_bams_list, ctrl_bams_list, majiq_gff3, threads, majiq_env, majiq_license, majiq_prob)
    return "MAJIQ pipeline run finished."

def execute_dexseq(pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, majiq_gff3, dexseq_env):
    dexseq_out = os.path.join(output_dir, "DEXSeq")
    if pipeline_mode == "Start from Raw FASTQ Files":
        exp_bams_list, ctrl_bams_list = get_bam_paths_from_fastq(exp_r1_input, ctrl_r1_input, output_dir)
    else:
        exp_bams_list, ctrl_bams_list = parse_bam_inputs(exp_bams, ctrl_bams)
        
    run_dexseq(dexseq_out, exp_bams_list, ctrl_bams_list, majiq_gff3, dexseq_env)
    return "DEXSeq pipeline run finished."

def execute_leafcutter(pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, leafcutter_env, threads):
    leafcutter_out = os.path.join(output_dir, "LeafCutter")
    if pipeline_mode == "Start from Raw FASTQ Files":
        exp_bams_list, ctrl_bams_list = get_bam_paths_from_fastq(exp_r1_input, ctrl_r1_input, output_dir)
    else:
        exp_bams_list, ctrl_bams_list = parse_bam_inputs(exp_bams, ctrl_bams)
        
    run_leafcutter(leafcutter_out, exp_bams_list, ctrl_bams_list, leafcutter_env, threads)
    return "LeafCutter pipeline run finished."

def execute_filtering(output_dir, pval_cutoff):
    rmats_out = os.path.join(output_dir, "rMATS")
    output_xlsx, total_events = process_rmats_output(rmats_out, pval_cutoff)
    if output_xlsx:
        return f"🎉 Extracted {total_events} significant events.", gr.update(value=output_xlsx, visible=True)
    else:
        return "No events matched your criteria or rMATS hasn't run yet.", gr.update(value=None, visible=False)

def check_majiq_download(output_dir):
    majiq_out = os.path.join(output_dir, "MAJIQ")
    majiq_filtered_tsv = os.path.join(majiq_out, "B1_vs_C_top_hits.tsv")
    if os.path.exists(majiq_filtered_tsv):
        return "MAJIQ TSV ready.", gr.update(value=majiq_filtered_tsv, visible=True)
    else:
        return "Run the MAJIQ module above to generate the TSV.", gr.update(value=None, visible=False)

# Toggle UI elements conditionally depending on mode selection
def toggle_mode_ui(mode):
    if mode == "Start from Raw FASTQ Files":
        return gr.update(visible=True), gr.update(visible=False), gr.update(visible=True), "## 3. Splicing Modules"
    else:
        return gr.update(visible=False), gr.update(visible=True), gr.update(visible=False), "## 2. Splicing Modules"

# ==========================================
# FRONTEND UI LAYOUT (Gradio Blocks)
# ==========================================
with gr.Blocks(title="Spliceome Pipeline") as demo:
    gr.Markdown("# 🧬 Modular RNA-Seq Splicing Pipeline")
    
    with gr.Sidebar():
        gr.Markdown("### Global Parameters")
        rmats_gtf = gr.Textbox(label="rMATS GTF", value="/home/alohani/majiq_scratch/gencode.v49.annotation.gtf")
        majiq_gff3 = gr.Textbox(label="MAJIQ GFF3", value="/home/alohani/majiq_scratch/gencode.v49.annotation.gff3")
        ref_fa = gr.Textbox(label="Reference Genome (.fa)", value="/home/data/hg38.fa")
        output_dir = gr.Textbox(label="Base Output Dir", value="/home/alohani/pipeline_output")
        threads = gr.Number(label="Threads", value=8, precision=0)
        read_length = gr.Number(label="Read Length", value=100, precision=0)
        
        gr.Markdown("### Filtering Thresholds")
        pval_cutoff = gr.Number(label="P-Value Cutoff (rMATS/DEXSeq)", value=0.05)
        majiq_prob = gr.Slider(label="MAJIQ Prob Cutoff", minimum=0.5, maximum=1.0, value=0.95, step=0.05)
        
        gr.Markdown("### Conda Environments")
        majiq_env = gr.Textbox(label="MAJIQ Env", value="majiq_scratch")
        dexseq_env = gr.Textbox(label="DEXSeq Env", value="/home/alohani/dexseq_env/")
        leafcutter_env = gr.Textbox(label="LeafCutter Env", value="leafcutter_env")
        edger_env = gr.Textbox(label="edgeR Env", value="edger_env")
        majiq_license = gr.Textbox(label="MAJIQ License", value="/home/alohani/majiq/license/majiq_license_academic_official.lic")

    # 1. INPUT FILES UI
    gr.Markdown("## 1. Pipeline Mode & Inputs")
    pipeline_mode = gr.Radio(
        choices=["Start from Raw FASTQ Files", "Start from Aligned BAM Files"], 
        value="Start from Raw FASTQ Files", 
        label="How are you starting your analysis?"
    )
    
    # Raw FASTQ Conditional View block
    with gr.Column(visible=True) as fastq_ui_block:
        gr.Info("Input your FASTQ files. They will be aligned and automatically sorted into edgeR/experimental and edgeR/control.")
        
        gr.Markdown("### Experimental Group")
        with gr.Row():
            exp_r1_input = gr.TextArea(label="Exp Read 1 (_R1)", lines=4)
            exp_r2_input = gr.TextArea(label="Exp Read 2 (_R2)", lines=4)
            
        gr.Markdown("### Control Group")
        with gr.Row():
            ctrl_r1_input = gr.TextArea(label="Ctrl Read 1 (_R1)", lines=4)
            ctrl_r2_input = gr.TextArea(label="Ctrl Read 2 (_R2)", lines=4)

    # BAM Paths Conditional View block
    with gr.Column(visible=False) as bam_ui_block:
        gr.Markdown("### BAM File Paths")
        with gr.Row():
            exp_bams = gr.TextArea(label="Experimental BAMs", lines=6)
            ctrl_bams = gr.TextArea(label="Control BAMs", lines=6)

    # 2. ALIGNMENT (Conditional UI section)
    with gr.Column(visible=True) as alignment_ui_block:
        gr.Markdown("## 2. Alignment Module")
        gr.Markdown("### Align & Count using Rsubread/edgeR")
        gen_index = gr.Checkbox(label="Generate Reference Index (.fa to index)", value=False)
        run_edger_btn = gr.Button("Run Full edgeR Pipeline", variant="primary")
        edger_status = gr.Textbox(label="edgeR Pipeline Log Status", interactive=False)
        
    # Dynamic header that updates based on the pipeline mode setting
    splicing_header = gr.Markdown("## 3. Splicing Modules")
    
    # 3. SPLICING MODULES (Tabs layout structure)
    with gr.Tabs():
        with gr.TabItem("rMATS"):
            run_rmats_btn = gr.Button("Run rMATS")
            rmats_status = gr.Textbox(label="rMATS Log Status", interactive=False)
            
        with gr.TabItem("MAJIQ"):
            run_majiq_btn = gr.Button("Run MAJIQ")
            majiq_status = gr.Textbox(label="MAJIQ Log Status", interactive=False)
            
        with gr.TabItem("DEXSeq"):
            run_dexseq_btn = gr.Button("Run DEXSeq Pipeline")
            dexseq_status = gr.Textbox(label="DEXSeq Log Status", interactive=False)
            
        with gr.TabItem("LeafCutter"):
            run_leafcutter_btn = gr.Button("Run LeafCutter")
            leafcutter_status = gr.Textbox(label="LeafCutter Log Status", interactive=False)

    # 4. FILTERING & DOWNLOADS
    gr.Markdown("## 4. Filter & Download Results")
    with gr.Row():
        with gr.Column():
            gr.Markdown("### rMATS Results")
            filter_rmats_btn = gr.Button("Filter & Generate rMATS Excel")
            rmats_filter_status = gr.Textbox(label="Filter Status", interactive=False)
            download_rmats_file = gr.File(label="📥 Download rMATS Excel", visible=False)
            
        with gr.Column():
            gr.Markdown("### MAJIQ Results")
            check_majiq_btn = gr.Button("Check & Fetch MAJIQ TSV")
            majiq_download_status = gr.Textbox(label="Download Status", interactive=False)
            download_majiq_file = gr.File(label="📥 Download MAJIQ TSV", visible=False)

    # ==========================================
    # EVENT TRIGGERS (Connecting Frontend to Backend)
    # ==========================================
    
    # Handle changing the pipeline selection mode dynamically
    pipeline_mode.change(
        fn=toggle_mode_ui,
        inputs=[pipeline_mode],
        outputs=[fastq_ui_block, bam_ui_block, alignment_ui_block, splicing_header]
    )
    
    # Trigger Pipeline Runs
    run_edger_btn.click(
        fn=execute_edger,
        inputs=[pipeline_mode, output_dir, exp_r1_input, exp_r2_input, ctrl_r1_input, ctrl_r2_input, ref_fa, threads, gen_index, edger_env],
        outputs=[edger_status]
    )
    run_rmats_btn.click(
        fn=execute_rmats,
        inputs=[pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, rmats_gtf, read_length, threads],
        outputs=[rmats_status]
    )
    run_majiq_btn.click(
        fn=execute_majiq,
        inputs=[pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, majiq_gff3, threads, majiq_env, majiq_license, majiq_prob],
        outputs=[majiq_status]
    )
    run_dexseq_btn.click(
        fn=execute_dexseq,
        inputs=[pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, majiq_gff3, dexseq_env],
        outputs=[dexseq_status]
    )
    run_leafcutter_btn.click(
        fn=execute_leafcutter,
        inputs=[pipeline_mode, output_dir, exp_r1_input, ctrl_r1_input, exp_bams, ctrl_bams, leafcutter_env, threads],
        outputs=[leafcutter_status]
    )
    
    # Trigger Downloader Systems
    filter_rmats_btn.click(
        fn=execute_filtering,
        inputs=[output_dir, pval_cutoff],
        outputs=[rmats_filter_status, download_rmats_file]
    )
    check_majiq_btn.click(
        fn=check_majiq_download,
        inputs=[output_dir],
        outputs=[majiq_download_status, download_majiq_file]
    )

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860)
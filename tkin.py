import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import threading
import sys
import os

# Import your existing logic
from edgeR_wrapper import run_edger_pipeline
from rmats import run_rmats
from majiq import run_majiq
from filter_results import process_rmats_output

class PipelineGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🧬 Modular RNA-Seq Splicing Pipeline (Native GUI)")
        self.root.geometry("900x700")

        # --- Create Scrollable Canvas ---
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Global Parameters
        ttk.Label(main_frame, text="Global Parameters", font=('Helvetica', 12, 'bold')).grid(row=0, column=0, sticky=tk.W)
        
        self.output_dir = self.create_input(main_frame, "Base Output Dir:", "/home/alohani/pipeline_output", 1)
        self.ref_fa = self.create_input(main_frame, "Reference Genome (.fa):", "/home/data/hg38.fa", 2)
        self.threads = self.create_input(main_frame, "Threads:", "22", 3)

        # FASTQ Inputs
        ttk.Label(main_frame, text="1. FASTQ Inputs", font=('Helvetica', 10, 'bold')).grid(row=4, column=0, pady=(10,0), sticky=tk.W)
        
        ttk.Label(main_frame, text="Experimental R1:").grid(row=5, column=0, sticky=tk.W)
        self.exp_r1_txt = scrolledtext.ScrolledText(main_frame, height=4, width=80)
        self.exp_r1_txt.grid(row=6, column=0, columnspan=2, pady=5)
        self.exp_r1_txt.insert(tk.END, "/home/alohani/QSIURP/MDAMB231/SRR20354100_1.fastq")

        ttk.Label(main_frame, text="Control R1:").grid(row=7, column=0, sticky=tk.W)
        self.ctrl_r1_txt = scrolledtext.ScrolledText(main_frame, height=4, width=80)
        self.ctrl_r1_txt.grid(row=8, column=0, columnspan=2, pady=5)
        self.ctrl_r1_txt.insert(tk.END, "/home/alohani/QSIURP/MDAMB231/SRR20354102_1.fastq")

        # Action Buttons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.grid(row=9, column=0, columnspan=2, pady=20)

        self.run_btn = ttk.Button(btn_frame, text="🚀 Run Alignment (edgeR)", command=self.start_edger_thread)
        self.run_btn.pack(side=tk.LEFT, padx=5)

        # Log Output
        ttk.Label(main_frame, text="Process Console Log:").grid(row=10, column=0, sticky=tk.W)
        self.log_area = scrolledtext.ScrolledText(main_frame, height=12, width=100, bg="black", fg="lime")
        self.log_area.grid(row=11, column=0, columnspan=2, pady=5)

    def create_input(self, parent, label, default, row):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky=tk.W, pady=2)
        var = tk.StringVar(value=default)
        entry = ttk.Entry(parent, textvariable=var, width=60)
        entry.grid(row=row, column=1, sticky=tk.W, pady=2)
        return var

    def log(self, message):
        self.log_area.insert(tk.END, message + "\n")
        self.log_area.see(tk.END)

    def start_edger_thread(self):
        # Disable button to prevent double clicks
        self.run_btn.config(state=tk.DISABLED)
        thread = threading.Thread(target=self.run_edger)
        thread.start()

    def run_edger(self):
        try:
            self.log("Starting pipeline... please wait.")
            # Gather inputs
            out_dir = os.path.join(self.output_dir.get(), "edgeR")
            ref = self.ref_fa.get()
            th = int(self.threads.get())
            
            # Parse text areas (simulating your r1/r2 logic)
            exp_r1 = [f.strip() for f in self.exp_r1_txt.get("1.0", tk.END).split('\n') if f.strip()]
            ctrl_r1 = [f.strip() for f in self.ctrl_r1_txt.get("1.0", tk.END).split('\n') if f.strip()]
            
            # Assuming R2 follows the _1 -> _2 naming convention
            exp_r2 = [f.replace("_1.fastq", "_2.fastq") for f in exp_r1]
            ctrl_r2 = [f.replace("_1.fastq", "_2.fastq") for f in ctrl_r1]

            # Run your actual wrapper
            # Note: gen_index is False here because we already built it!
            run_edger_pipeline(out_dir, exp_r1, exp_r2, ctrl_r1, ctrl_r2, ref, th, False, "edger_env")
            
            self.log("✅ edgeR Pipeline Completed Successfully!")
            messagebox.showinfo("Success", "Alignment and Counting finished!")
        except Exception as e:
            self.log(f"❌ ERROR: {str(e)}")
            messagebox.showerror("Error", str(e))
        finally:
            self.run_btn.config(state=tk.NORMAL)

if __name__ == "__main__":
    root = tk.Tk()
    gui = PipelineGUI(root)
    root.mainloop()
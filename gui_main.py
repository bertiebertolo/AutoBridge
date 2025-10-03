#!/usr/bin/env python3
"""
AutoBLAST-GUI: Main GUI entrypoint (Tkinter)

Three tabs:
 - Ortholog mapping (Ensembl + RBH fallback)
 - Protein BLAST (blastp)
 - Gene BLAST (blastn / tblastn style)

This is a simple, extendable GUI that calls functions in ensembl_rbh.py and blast_utils.py.
"""
import os
import tkinter as tk
from tkinter import ttk, filedialog, simpledialog, messagebox
import threading

from ensembl_rbh import map_with_ensembl_and_rbh_ncbi, map_genes  # import primary functions
from blast_utils import check_blast_tools, run_blastp, run_blastn, make_blast_db

APP_TITLE = "AutoBLAST-GUI"

class App:
    def __init__(self, root):
        self.root = root
        root.title(APP_TITLE)
        self.nb = ttk.Notebook(root)
        self.nb.pack(fill="both", expand=True)

        self._build_ortholog_tab()
        self._build_protein_blast_tab()
        self._build_gene_blast_tab()

    def _build_ortholog_tab(self):
        frame = ttk.Frame(self.nb)
        self.nb.add(frame, text="Ortholog Mapping")

        lbl = ttk.Label(frame, text="Input CSV (GENE_ID, HD, JO):")
        lbl.grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.orth_input_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.orth_input_var, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(frame, text="Browse", command=self._pick_orth_input).grid(row=0, column=2, padx=5)

        ttk.Label(frame, text="Output CSV:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.orth_output_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.orth_output_var, width=60).grid(row=1, column=1, padx=5)
        ttk.Button(frame, text="Choose", command=self._pick_orth_output).grid(row=1, column=2, padx=5)

        ttk.Label(frame, text="NCBI Auto-download proteomes?").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.orth_auto_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frame, variable=self.orth_auto_var).grid(row=2, column=1, sticky="w")

        ttk.Label(frame, text="Ensembl target species substring:").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.orth_target_var = tk.StringVar(value="culex")
        ttk.Entry(frame, textvariable=self.orth_target_var, width=30).grid(row=3, column=1, sticky="w", padx=5)

        ttk.Label(frame, text="Entrez email (for NCBI):").grid(row=4, column=0, sticky="w", padx=5, pady=5)
        self.orth_email_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.orth_email_var, width=40).grid(row=4, column=1, padx=5)

        self.orth_progress = tk.StringVar(value="")
        ttk.Label(frame, textvariable=self.orth_progress).grid(row=5, column=0, columnspan=3, padx=5, pady=5)

        ttk.Button(frame, text="Run Mapping", command=self._run_mapping_thread).grid(row=6, column=1, pady=8)

    def _build_protein_blast_tab(self):
        frame = ttk.Frame(self.nb)
        self.nb.add(frame, text="Protein BLAST (blastp)")

        ttk.Label(frame, text="Query protein FASTA:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.pblast_query = tk.StringVar()
        ttk.Entry(frame, textvariable=self.pblast_query, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(frame, text="Browse", command=lambda: self._pick_file(self.pblast_query, "Select query FASTA")).grid(row=0, column=2)

        ttk.Label(frame, text="Subject DB FASTA (makeblastdb will be used):").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.pblast_db = tk.StringVar()
        ttk.Entry(frame, textvariable=self.pblast_db, width=60).grid(row=1, column=1, padx=5)
        ttk.Button(frame, text="Browse", command=lambda: self._pick_file(self.pblast_db, "Select subject FASTA")).grid(row=1, column=2)

        ttk.Label(frame, text="Output Path:").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.pblast_out = tk.StringVar()
        ttk.Entry(frame, textvariable=self.pblast_out, width=60).grid(row=2, column=1, padx=5)
        ttk.Button(frame, text="Choose", command=lambda: self._pick_file(self.pblast_out, "Save output", save=True)).grid(row=2, column=2)

        ttk.Button(frame, text="Run blastp", command=self._run_blastp_thread).grid(row=3, column=1, pady=10)
        self.pblast_status = tk.StringVar(value="")
        ttk.Label(frame, textvariable=self.pblast_status).grid(row=4, column=0, columnspan=3, padx=5)

    def _build_gene_blast_tab(self):
        frame = ttk.Frame(self.nb)
        self.nb.add(frame, text="Gene BLAST (blastn / tblastn)")

        ttk.Label(frame, text="Query nucleotide FASTA:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.gblast_query = tk.StringVar()
        ttk.Entry(frame, textvariable=self.gblast_query, width=60).grid(row=0, column=1, padx=5)
        ttk.Button(frame, text="Browse", command=lambda: self._pick_file(self.gblast_query, "Select query FASTA")).grid(row=0, column=2)

        ttk.Label(frame, text="Subject DB FASTA (makeblastdb -dbtype nucl):").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.gblast_db = tk.StringVar()
        ttk.Entry(frame, textvariable=self.gblast_db, width=60).grid(row=1, column=1, padx=5)
        ttk.Button(frame, text="Browse", command=lambda: self._pick_file(self.gblast_db, "Select subject FASTA")).grid(row=1, column=2)

        ttk.Label(frame, text="Output Path:").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.gblast_out = tk.StringVar()
        ttk.Entry(frame, textvariable=self.gblast_out, width=60).grid(row=2, column=1, padx=5)
        ttk.Button(frame, text="Choose", command=lambda: self._pick_file(self.gblast_out, "Save output", save=True)).grid(row=2, column=2)

        ttk.Button(frame, text="Run gene BLAST", command=self._run_blastn_thread).grid(row=3, column=1, pady=10)
        self.gblast_status = tk.StringVar(value="")
        ttk.Label(frame, textvariable=self.gblast_status).grid(row=4, column=0, columnspan=3, padx=5)

    # File pickers
    def _pick_orth_input(self):
        path = filedialog.askopenfilename(title="Select input CSV", filetypes=[("CSV","*.csv"),("All files","*.*")])
        if path:
            self.orth_input_var.set(path)

    def _pick_orth_output(self):
        path = filedialog.asksaveasfilename(title="Save output CSV", defaultextension=".csv", filetypes=[("CSV","*.csv")])
        if path:
            self.orth_output_var.set(path)

    def _pick_file(self, var, title, save=False):
        if save:
            p = filedialog.asksaveasfilename(title=title, defaultextension=".txt")
        else:
            p = filedialog.askopenfilename(title=title)
        if p:
            var.set(p)

    # Threaded runners (so GUI stays responsive)
    def _run_mapping_thread(self):
        infile = self.orth_input_var.get().strip()
        outfile = self.orth_output_var.get().strip()
        email = self.orth_email_var.get().strip()
        auto = self.orth_auto_var.get()
        target = self.orth_target_var.get().strip()
        if not infile or not outfile or not email:
            messagebox.showwarning("Missing", "Please select input, output and provide an email.")
            return
        def worker():
            try:
                self.orth_progress.set("Running mapping...")
                out, n = map_with_ensembl_and_rbh_ncbi(infile, outfile, email, target, auto, sleep=0.25)
                self.orth_progress.set(f"Completed: wrote {n} rows to {out}")
            except Exception as e:
                self.orth_progress.set(f"Error: {e}")
        threading.Thread(target=worker, daemon=True).start()

    def _run_blastp_thread(self):
        q = self.pblast_query.get().strip()
        db = self.pblast_db.get().strip()
        out = self.pblast_out.get().strip()
        if not q or not db or not out:
            messagebox.showwarning("Missing", "Please choose query, subject FASTA and output")
            return
        def worker():
            try:
                self.pblast_status.set("Checking BLAST tools...")
                check_blast_tools()
                self.pblast_status.set("Building DB...")
                db_base = os.path.splitext(db)[0]
                make_blast_db(db, db_base)
                self.pblast_status.set("Running blastp...")
                run_blastp(q, db_base, out)
                self.pblast_status.set(f"blastp complete -> {out}")
            except Exception as e:
                self.pblast_status.set(f"Error: {e}")
        threading.Thread(target=worker, daemon=True).start()

    def _run_blastn_thread(self):
        q = self.gblast_query.get().strip()
        db = self.gblast_db.get().strip()
        out = self.gblast_out.get().strip()
        if not q or not db or not out:
            messagebox.showwarning("Missing", "Please choose query, subject FASTA and output")
            return
        def worker():
            try:
                self.gblast_status.set("Checking BLAST tools...")
                check_blast_tools()
                self.gblast_status.set("Building DB (nucl)...")
                db_base = os.path.splitext(db)[0]
                # makeblastdb for nucleotide
                subprocess_cmd = f"makeblastdb -in \"{db}\" -dbtype nucl -out \"{db_base}\""
                os.system(subprocess_cmd)
                self.gblast_status.set("Running blastn...")
                run_blastn(q, db_base, out)
                self.gblast_status.set(f"blastn complete -> {out}")
            except Exception as e:
                self.gblast_status.set(f"Error: {e}")
        threading.Thread(target=worker, daemon=True).start()

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
#!/usr/bin/env python3
"""
AUTOBridge GUI — NCBI-first suggest + mapping (no local FASTA inputs)

This GUI uses NCBI Taxonomy for species suggestion and supports an NCBI-only RBH
mapping flow (downloads proteomes via ncbi_utils and runs BLAST-based RBH).
"""
import os
import threading
import tkinter as tk
from tkinter import ttk, filedialog, simpledialog, messagebox
import requests
import time
import subprocess

from ensembl_rbh import map_with_ensembl_and_rbh_ncbi, map_genes_ncbi_rbh  # core functions
from blast_utils import check_blast_tools, run_blastp, run_blastn, make_blast_db

APP_TITLE = "AUTOBridge"
ENSEMBL_REST = "https://rest.ensembl.org"
DEFAULT_LOG = "autobridge.log"


class App:
    def __init__(self, root):
        self.root = root
        root.title(APP_TITLE)
        self.nb = ttk.Notebook(root)
        self.nb.pack(fill="both", expand=True)

        self._build_ortholog_tab()
        self._build_protein_blast_tab()
        self._build_gene_blast_tab()

        # log tail control
        self._log_tail_thread = None
        self._log_tail_stop = threading.Event()

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

        # Source/Target taxa
        ttk.Label(frame, text="Source taxon (e.g. 'Aedes aegypti'):", ).grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.source_taxon_var = tk.StringVar(value="Aedes aegypti")
        ttk.Entry(frame, textvariable=self.source_taxon_var, width=40).grid(row=2, column=1, padx=5)

        ttk.Label(frame, text="Target taxon (e.g. 'Culex pipiens pallens'):", wraplength=400).grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.orth_target_var = tk.StringVar(value="Culex pipiens pallens")
        ttk.Entry(frame, textvariable=self.orth_target_var, width=40).grid(row=3, column=1, padx=5)
        ttk.Button(frame, text="Suggest species (NCBI)", command=self._suggest_species).grid(row=3, column=2, padx=5)

        # Auto-download and mode
        self.orth_auto_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frame, variable=self.orth_auto_var, text="Auto-download from NCBI").grid(row=4, column=1, sticky="w", padx=5, pady=3)
        self.use_ncbi_rbh_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frame, variable=self.use_ncbi_rbh_var, text="Use NCBI RBH (skip Ensembl)").grid(row=4, column=2, sticky="w", padx=5, pady=3)

        # Entrez email
        ttk.Label(frame, text="Entrez email (for NCBI):").grid(row=5, column=0, sticky="w", padx=5, pady=5)
        self.orth_email_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.orth_email_var, width=40).grid(row=5, column=1, padx=5)

        guidance = ("This mode downloads proteomes from NCBI datasets and runs BLAST-based RBH to infer orthologs.\n"
                    "Ensure you have BLAST+ installed and ncbi-datasets-cli available for best results.")
        ttk.Label(frame, text=guidance, foreground="gray", wraplength=700).grid(row=6, column=0, columnspan=3, padx=5, pady=6)

        # Show progress checkbox
        self.show_progress_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frame, text="Show progress output", variable=self.show_progress_var).grid(row=7, column=1, sticky="w", padx=5, pady=3)

        # Single-line status + log viewer
        self.orth_progress = tk.StringVar(value="")
        ttk.Label(frame, textvariable=self.orth_progress).grid(row=8, column=0, columnspan=3, padx=5, pady=5)

        ttk.Label(frame, text="Live log (autobridge.log):").grid(row=9, column=0, sticky="nw", padx=5)
        self.log_text = tk.Text(frame, width=90, height=12, wrap="none")
        self.log_text.grid(row=9, column=1, columnspan=2, padx=5, pady=5)
        self.log_text.configure(state="disabled")
        log_sb = ttk.Scrollbar(frame, orient="vertical", command=self._scroll_log)
        log_sb.grid(row=9, column=3, sticky="ns")
        self.log_text['yscrollcommand'] = log_sb.set

        ttk.Button(frame, text="Run Mapping", command=self._run_mapping_thread).grid(row=10, column=1, pady=8)

    # protein/gene blast tabs are unchanged
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

        self.pblast_status = tk.StringVar(value="")
        ttk.Button(frame, text="Run blastp", command=self._run_blastp_thread).grid(row=3, column=1, pady=10)
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

        self.gblast_status = tk.StringVar(value="")
        ttk.Button(frame, text="Run gene BLAST", command=self._run_blastn_thread).grid(row=3, column=1, pady=10)
        ttk.Label(frame, textvariable=self.gblast_status).grid(row=4, column=0, columnspan=3, padx=5)

    # File pickers
    def _pick_orth_input(self):
        path = filedialog.askopenfilename(title="Select input CSV", filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if path:
            self.orth_input_var.set(path)

    def _pick_orth_output(self):
        path = filedialog.asksaveasfilename(title="Save output CSV", defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if path:
            self.orth_output_var.set(path)

    def _pick_file(self, var, title, save=False):
        if save:
            p = filedialog.asksaveasfilename(title=title, defaultextension=".txt")
        else:
            p = filedialog.askopenfilename(title=title)
        if p:
            var.set(p)

    # Suggest species using NCBI Taxonomy (esearch + esummary)
    def _suggest_species(self):
        import difflib
        sub = simpledialog.askstring("Suggest species", "Enter a substring to search NCBI taxonomy (e.g. 'culex', 'pipiens', 'human'):")
        if not sub:
            return

        email = self.orth_email_var.get().strip()
        try:
            esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {"db": "taxonomy", "term": sub, "retmax": 200, "retmode": "json"}
            if email:
                params["email"] = email
            r = requests.get(esearch_url, params=params, timeout=30)
            r.raise_for_status()
            j = r.json()
            idlist = j.get("esearchresult", {}).get("idlist", [])
            matches = []
            if idlist:
                esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                params2 = {"db": "taxonomy", "id": ",".join(idlist), "retmode": "json"}
                if email:
                    params2["email"] = email
                r2 = requests.get(esummary_url, params=params2, timeout=30)
                r2.raise_for_status()
                summary = r2.json().get("result", {})
                uids = summary.get("uids", [])
                for uid in uids:
                    item = summary.get(uid, {})
                    sci = item.get("scientificname", "") or ""
                    common = item.get("commonname", "") or ""
                    rank = item.get("rank", "") or ""
                    disp = sci
                    if common:
                        disp = f"{sci} ({common})"
                    if rank:
                        disp = f"{disp} — {rank}"
                    matches.append((sci, disp))
            if matches:
                # dedupe while preserving order
                seen = set()
                uniq = []
                for sci, disp in matches:
                    if sci not in seen:
                        seen.add(sci)
                        uniq.append((sci, disp))
                self._show_species_selection(uniq, title="NCBI Taxonomy matches (double-click to select)")
                return
        except Exception:
            pass

        # fallback to Ensembl names (rare; kept for resilience)
        try:
            url = ENSEMBL_REST + "/info/species"
            headers = {"Content-Type": "application/json", "User-Agent": APP_TITLE}
            r = requests.get(url, headers=headers, timeout=30)
            r.raise_for_status()
            data = r.json()
            matches = []
            for s in data.get("species", []):
                name = s.get("name", "")
                disp = s.get("display_name", "")
                if sub.lower() in name.lower() or sub.lower() in disp.lower():
                    matches.append((name, disp))
            if not matches:
                names = [s.get("name", "") for s in data.get("species", []) if s.get("name", "")]
                close = difflib.get_close_matches(sub.lower(), [n.lower() for n in names], n=10, cutoff=0.4)
                if close:
                    candidates = [(n, "") for n in close]
                    self._show_species_selection(candidates, title="Close matches - pick one")
                else:
                    messagebox.showinfo("No matches", f"No NCBI or Ensembl species matched '{sub}'. Try a different substring.")
            else:
                self._show_species_selection(matches, title="Ensembl species matches (double-click to select)")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to fetch species list from NCBI and Ensembl: {e}")

    def _show_species_selection(self, name_disp_pairs, title="Select species"):
        win = tk.Toplevel(self.root)
        win.title(title)
        listbox = tk.Listbox(win, width=80, height=20)
        for name, disp in name_disp_pairs:
            text = f"{name} — {disp}" if disp else name
            listbox.insert("end", text)
        listbox.pack(side="left", fill="both", expand=True)
        sb = ttk.Scrollbar(win, command=listbox.yview)
        sb.pack(side="right", fill="y")
        listbox.config(yscrollcommand=sb.set)

        def on_select(event=None):
            sel = listbox.curselection()
            if not sel:
                return
            text = listbox.get(sel[0])
            canonical = text.split(" — ", 1)[0]
            # populate both source and target taxon with selected scientific name if empty
            if not self.source_taxon_var.get().strip():
                self.source_taxon_var.set(canonical)
            else:
                self.orth_target_var.set(canonical)
            win.destroy()

        listbox.bind("<Double-Button-1>", on_select)
        ttk.Button(win, text="Select", command=on_select).pack(pady=6)

    def _scroll_log(self, *args):
        self.log_text.yview(*args)

    def _start_log_tail(self, log_path: str):
        """Start a background thread to tail log_path and write to the log_text widget."""
        # stop any existing tail
        self._stop_log_tail()

        self._log_tail_stop.clear()

        def tail():
            last_size = 0
            try:
                while not self._log_tail_stop.is_set():
                    if os.path.exists(log_path):
                        try:
                            size = os.path.getsize(log_path)
                            if size < last_size:
                                # file rotated/truncated
                                last_size = 0
                            if size > last_size:
                                with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
                                    fh.seek(last_size)
                                    new = fh.read()
                                last_size = os.path.getsize(log_path)
                                if new:
                                    self._append_log_text(new)
                        except Exception:
                            pass
                    time.sleep(1.0)
            except Exception:
                pass

        self._log_tail_thread = threading.Thread(target=tail, daemon=True)
        self._log_tail_thread.start()

    def _stop_log_tail(self):
        if self._log_tail_thread and self._log_tail_thread.is_alive():
            self._log_tail_stop.set()
            self._log_tail_thread.join(timeout=2.0)
        self._log_tail_thread = None
        self._log_tail_stop.clear()

    def _append_log_text(self, text: str):
        def append():
            self.log_text.configure(state="normal")
            self.log_text.insert("end", text)
            self.log_text.see("end")
            self.log_text.configure(state="disabled")
        self.root.after(0, append)

    def _run_mapping_thread(self):
        infile = self.orth_input_var.get().strip()
        outfile = self.orth_output_var.get().strip()
        email = self.orth_email_var.get().strip()
        auto = self.orth_auto_var.get()
        target = self.orth_target_var.get().strip()
        source = self.source_taxon_var.get().strip()
        show_progress = self.show_progress_var.get()
        use_ncbi_rbh = self.use_ncbi_rbh_var.get()

        if not infile or not outfile or not email:
            messagebox.showwarning("Missing", "Please select input, output and provide an email.")
            return

        # clear log viewer
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.configure(state="disabled")

        def worker():
            try:
                # start tailing the log
                self.orth_progress.set("Starting mapping run...")
                self._start_log_tail(DEFAULT_LOG)

                if use_ncbi_rbh:
                    self.orth_progress.set(f"Running NCBI RBH mapping: {source} -> {target}")
                    out, n = map_genes_ncbi_rbh(infile, outfile, email, source, target,
                                                sleep=0.25, show_progress=show_progress, log_file=DEFAULT_LOG)
                else:
                    # default: call Ensembl+RBH wrapper (keeps backward compatibility)
                    out, n = map_with_ensembl_and_rbh_ncbi(infile, outfile, email, target,
                                                          auto, sleep=0.25, show_progress=show_progress, log_file=DEFAULT_LOG)
                self.orth_progress.set(f"Completed: wrote {n} rows to {out}")
            except Exception as e:
                self.orth_progress.set(f"Error: {e}")
                self._append_log_text(f"\nERROR: {e}\n")
            finally:
                self._stop_log_tail()
        threading.Thread(target=worker, daemon=True).start()

    # Blast helper threads (unchanged)
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
                try:
                    make_blast_db(db, db_base)
                except TypeError:
                    subprocess_cmd = f'makeblastdb -in "{db}" -dbtype prot -out "{db_base}"'
                    os.system(subprocess_cmd)
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
                try:
                    make_blast_db(db, db_base, dbtype="nucl")
                except TypeError:
                    subprocess_cmd = f'makeblastdb -in "{db}" -dbtype nucl -out "{db_base}"'
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

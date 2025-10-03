#!/usr/bin/env python3
"""
Map Aedes aegypti AAEL gene IDs to Culex pipiens pallens orthologs using Ensembl REST
(primary) plus reciprocal-BLAST (RBH) fallback for high accuracy.

Features
- GUI file pickers to choose input CSV and output file.
- Uses Ensembl/VectorBase homology endpoint first (accepts one-to-one orthologs).
- For missing/ambiguous cases, performs reciprocal best-hit (RBH) using BLAST+:
  - Requires local Culex proteome FASTA and (preferably) local Aedes proteome FASTA.
  - If a local Aedes proteome isn't supplied, the script can fetch each AAEL protein
    sequence from Ensembl REST (slower).
- Outputs CSV with details: query_id, target_gene_id, method_used, orthology_type,
  percent_identity, coverage, evalue, HD, JO, note.

Requirements
- BLAST+ (makeblastdb, blastp) installed and on PATH.
- Python packages: pandas, requests, biopython
    pip install pandas requests biopython

Recommended proteome sources (you may download and supply local FASTA files):
- Aedes aegypti (VectorBase / Ensembl Genomes): ftp://ftp.ensemblgenomes.org/pub/metazoa/
- Culex pipiens (VectorBase / NCBI): check VectorBase/NCBI FTPs for the species. Filename
  usually contains "pep" or "protein" for proteomes.

Usage
- Run: python ensembl_map_aael_to_culex_rbh.py
- Follow dialogs (select input CSV, choose output path, provide email).
"""

import os
import sys
import csv
import time
import subprocess
import tempfile
import shutil
import requests
import pandas as pd
from io import StringIO
from collections import defaultdict
from typing import Dict, Optional, Tuple, List
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ==== Configuration defaults ====
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS_JSON = {"Accept": "application/json"}
HEADERS_FASTA = {"Accept": "text/x-fasta"}
BLAST_EVALUE = 1e-5
BLAST_IDENTITY_MIN = 30.0      # percent
BLAST_COVERAGE_MIN = 0.5       # fraction of query covered
BLAST_MAX_TARGETS = 5

# ==== Utilities ====
def which(exe_name: str) -> Optional[str]:
    return shutil.which(exe_name)

def check_blast_tools():
    for exe in ("makeblastdb", "blastp"):
        if which(exe) is None:
            raise FileNotFoundError(f"{exe} not found on PATH. Please install BLAST+ and ensure {exe} is on PATH.")

def make_blast_db(fasta_path: str, db_base: str):
    """Create BLAST protein DB at db_base (skip if files exist)."""
    # check for db files (blastdb will create .pin/.phr/.psq or .p* depending on version)
    if os.path.exists(db_base + ".pin") or os.path.exists(db_base + ".psq") or os.path.exists(db_base + ".phr"):
        return db_base
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "prot", "-out", db_base]
    subprocess.run(cmd, check=True)
    return db_base

def run_blastp_query(query_fasta: str, db_base: str, outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
    cmd = [
        "blastp", "-query", query_fasta, "-db", db_base,
        "-outfmt", outfmt,
        "-max_target_seqs", str(BLAST_MAX_TARGETS),
        "-evalue", str(BLAST_EVALUE)
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, text=True)
    if proc.returncode not in (0, 1):  # blastp returns 1 if no hits found
        raise RuntimeError(f"blastp failed: {proc.stderr.strip()}")
    return proc.stdout

def parse_blast_tab6_line(line: str):
    fields = line.strip().split('\t')
    # returns qid, sid, pident, align_len, evalue
    return {
        "qseqid": fields[0],
        "sseqid": fields[1],
        "pident": float(fields[2]),
        "length": int(fields[3]),
        "evalue": float(fields[10]),
    }

# ==== Ensembl REST helpers ====
def ensembl_homology(gene_id: str, session: requests.Session, user_agent: str, target_species_substr: str):
    url = f"{ENSEMBL_REST}/homology/id/{gene_id}"
    params = {"type": "orthologues", "sequence": 0, "format": "full"}
    headers = HEADERS_JSON.copy()
    headers["User-Agent"] = user_agent
    resp = session.get(url, params=params, headers=headers, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    rows = []
    for rec in data.get("data", []):
        for h in rec.get("homologies", []):
            tgt = h.get("target", {})
            species = tgt.get("species", "")
            if target_species_substr.lower() not in species.lower():
                continue
            target_gene = tgt.get("gene_id") or tgt.get("id") or ""
            rows.append({
                "target_gene_id": target_gene,
                "target_protein_id": tgt.get("id", ""),
                "target_species": species,
                "orthology_type": h.get("type", ""),
                "confidence": h.get("confidence", ""),
                "perc_id": h.get("percentage_identities") or h.get("perc_id") or ""
            })
    return rows

def ensembl_fetch_protein_sequence(identifier: str, session: requests.Session, user_agent: str) -> Optional[SeqRecord]:
    """Fetch protein sequence (FASTA) for the given stable identifier via Ensembl REST.
    Returns SeqRecord or None on failure.
    """
    # Try protein sequence by id
    url = f"{ENSEMBL_REST}/sequence/id/{identifier}"
    params = {"type": "protein"}
    headers = HEADERS_FASTA.copy()
    headers["User-Agent"] = user_agent
    resp = session.get(url, params=params, headers=headers, timeout=30)
    if resp.status_code == 200:
        fasta = resp.text
        try:
            rec = next(SeqIO.parse(StringIO(fasta), "fasta"))
            return rec
        except Exception:
            return None
    else:
        return None

# ==== FASTA utilities ====
def index_fasta_by_id_substring(fasta_path: str) -> Dict[str, SeqRecord]:
    """Index fasta records by header and also keep a mapping from gene-id-like substrings to records.
    We'll check if query gene id (AAEL...) is substring of record.id or record.description.
    """
    index = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        index[rec.id] = rec
    return index

def find_records_containing_substr(index: Dict[str, SeqRecord], substr: str) -> List[SeqRecord]:
    substr = substr.strip()
    out = []
    for key, rec in index.items():
        if substr in key or substr in rec.description:
            out.append(rec)
    return out

# ==== RBH logic ====
def run_rbh_for_query(query_gene: str, aedes_index: Optional[Dict[str, SeqRecord]],
                      culex_db_base: str, aedes_db_base: Optional[str],
                      session: requests.Session, user_agent: str,
                      tmpdir: str) -> Tuple[Optional[str], Dict]:
    """
    Returns (target_culex_gene_id_or_protein_id, metrics dict)
    metrics includes method_used='rbh' or 'rbh_no_hit' or 'rbh_ambiguous', percent_identity, coverage, evalue
    """
    # 1) obtain query protein sequence (prefer local index)
    seqrec = None
    if aedes_index is not None:
        matches = find_records_containing_substr(aedes_index, query_gene)
        if matches:
            # prefer exact matching record id containing gene string; otherwise take first
            seqrec = matches[0]
    if seqrec is None:
        # fallback to Ensembl REST for this AAEL id
        seqrec = ensembl_fetch_protein_sequence(query_gene, session, user_agent)
    if seqrec is None:
        return None, {"method_used": "rbh_no_query_seq", "note": "no_query_sequence"}

    # write query to temp fasta
    qf = os.path.join(tmpdir, f"query_{query_gene}.fa")
    SeqIO.write(seqrec, qf, "fasta")

    # BLASTP query -> culex DB
    out = run_blastp_query(qf, culex_db_base)
    if not out.strip():
        return None, {"method_used": "rbh_no_hit", "note": "no_hit_to_culex"}

    # parse top hit
    top_line = out.strip().splitlines()[0]
    top = top_line.split('\t')
    sseqid = top[1]
    pident = float(top[2])
    align_len = int(top[3])
    evalue = float(top[10])
    # compute coverage
    qlen = len(seqrec.seq)
    coverage = align_len / qlen if qlen > 0 else 0.0

    # quality filters
    if pident < BLAST_IDENTITY_MIN or coverage < BLAST_COVERAGE_MIN or evalue > BLAST_EVALUE:
        return None, {"method_used": "rbh_hit_low_quality", "pident": pident, "coverage": coverage, "evalue": evalue, "note": "hit_below_thresholds"}

    # get the subject sequence from Culex fasta (we need the sequence to BLAST back).
    # We assume culex_db_base + ".fa" is not present; instead we need the original culex FASTA path.
    # To avoid needing to re-open DB internals, we will require that culex_db_base was built from a provided fasta
    # and that we stored that path in culex_db_base + ".source_fasta" (we create this earlier).
    source_info_path = culex_db_base + ".source_fasta"
    if not os.path.exists(source_info_path):
        return None, {"method_used": "rbh_no_culex_fasta_info", "note": "culex_source_not_recorded"}
    with open(source_info_path) as fh:
        culex_fasta = fh.read().strip()
    # find record
    culex_seqrec = None
    for rec in SeqIO.parse(culex_fasta, "fasta"):
        if sseqid == rec.id or sseqid in rec.id or sseqid in rec.description:
            culex_seqrec = rec
            break
    if culex_seqrec is None:
        # fallback: use blastdbcmd to fetch sequence (if available)
        if which("blastdbcmd"):
            try:
                cmd = ["blastdbcmd", "-db", culex_db_base, "-entry", sseqid, "-outfmt", "%f"]
                res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
                culex_seqrec = next(SeqIO.parse(StringIO(res.stdout), "fasta"))
            except Exception:
                return None, {"method_used": "rbh_no_subject_seq", "note": "cannot_retrieve_subject_seq"}
        else:
            return None, {"method_used": "rbh_no_subject_seq_no_blastdbcmd", "note": "cannot_retrieve_subject_seq"}

    # write subject sequence to temp fasta
    subj_fa = os.path.join(tmpdir, f"subj_{sseqid}.fa")
    SeqIO.write(culex_seqrec, subj_fa, "fasta")

    # BLAST subject -> aedes DB
    if aedes_db_base is None:
        # we can BLAST subject against local aedes fasta by making a DB on the fly if we have aedes_index source
        if aedes_index is None:
            return None, {"method_used": "rbh_no_aedes_db", "note": "no_aedes_db_available_for_reciprocal"}
        # write aedes index to temporary fasta DB
        aedes_tmp_fa = os.path.join(tmpdir, "aedes_from_index.fa")
        SeqIO.write(list(aedes_index.values()), aedes_tmp_fa, "fasta")
        aedes_db_base_tmp = os.path.join(tmpdir, "aedes_tmp_db")
        make_blast_db(aedes_tmp_fa, aedes_db_base_tmp)
        aedes_db_for_rbh = aedes_db_base_tmp
    else:
        aedes_db_for_rbh = aedes_db_base

    out2 = run_blastp_query(subj_fa, aedes_db_for_rbh)
    if not out2.strip():
        return None, {"method_used": "rbh_no_reciprocal_hit", "note": "no_reciprocal_hit"}
    top2 = out2.strip().splitlines()[0].split('\t')
    top2_qid = top2[0]
    top2_sseqid = top2[1]
    # For reciprocal best hit, we want the top hit's subject (subject from culex->aedes blast)
    # to map back to the original AAEL identifier. Depending on DB contents, top2_sseqid or top2_qid might be used.
    # We check if original query_gene appears in top2_sseqid or top2_qid or in the description of the aedes record.
    # We'll read the top hit record's id and check substring.
    reciprocal_hit_id = top2[1]
    # Check if reciprocal_hit_id or top2_qid contains query_gene
    if query_gene in reciprocal_hit_id or query_gene in top2_qid:
        return reciprocal_hit_id, {"method_used": "rbh", "pident": pident, "coverage": coverage, "evalue": evalue}
    else:
        return None, {"method_used": "rbh_nonreciprocal", "pident": pident, "coverage": coverage, "evalue": evalue, "reciprocal_top": reciprocal_hit_id}

# ==== Main mapping flow ====
def map_with_ensembl_and_rbh(input_csv: str, output_csv: str, email: str,
                             target_species: str, culex_fasta_path: Optional[str],
                             aedes_fasta_path: Optional[str],
                             sleep: float = 0.25):
    check_blast_tools()
    session = requests.Session()
    user_agent = f"ensembl-rbh-script ({email})" if email else "ensembl-rbh-script"

    # load input CSV
    df = pd.read_csv(input_csv)
    if "GENE_ID" not in df.columns:
        raise ValueError("Input CSV must contain GENE_ID column")
    # prepare BLAST DBs if FASTAs provided
    tmpdir_global = tempfile.mkdtemp(prefix="ensembl_rbh_")
    try:
        culex_db_base = None
        if culex_fasta_path:
            culex_db_base = os.path.join(tmpdir_global, "culex_db")
            make_blast_db(culex_fasta_path, culex_db_base)
            # record source fasta path for later retrieval
            with open(culex_db_base + ".source_fasta", "w") as fh:
                fh.write(os.path.abspath(culex_fasta_path))
        aedes_db_base = None
        aedes_index = None
        if aedes_fasta_path:
            aedes_db_base = os.path.join(tmpdir_global, "aedes_db")
            make_blast_db(aedes_fasta_path, aedes_db_base)
            aedes_index = index_fasta_by_id_substring(aedes_fasta_path)

        # results
        out_rows = []
        total = len(df)
        for i, row in df.iterrows():
            qid = str(row["GENE_ID"]).strip().strip('"')
            hd = row.get("HD", "")
            jo = row.get("JO", "")
            # 1) try Ensembl homology
            try:
                homs = ensembl_homology(qid, session, user_agent, target_species)
            except Exception as e:
                homs = []
                ensembl_error = str(e)
            else:
                ensembl_error = ""

            chosen_rows = []
            if homs:
                # prefer one-to-one
                for h in homs:
                    typ = h.get("orthology_type", "").lower()
                    if "one2one" in typ or "one_to_one" in typ or "one-to-one" in typ:
                        chosen_rows.append((h, "ensembl_one2one"))
                if not chosen_rows:
                    # accept others
                    for h in homs:
                        chosen_rows.append((h, "ensembl_other"))

            if chosen_rows:
                for h, method_tag in chosen_rows:
                    out_rows.append({
                        "query_id": qid,
                        "target_gene_id": h.get("target_gene_id", ""),
                        "method_used": "ensembl",
                        "orthology_type": h.get("orthology_type", ""),
                        "percent_identity": h.get("perc_id", ""),
                        "coverage": "",
                        "evalue": "",
                        "HD": hd,
                        "JO": jo,
                        "note": method_tag
                    })
            else:
                # need RBH fallback for accuracy
                if culex_db_base is None:
                    # we cannot do RBH without Culex DB; record Ensembl attempt and mark no_culex_db
                    out_rows.append({
                        "query_id": qid,
                        "target_gene_id": "",
                        "method_used": "ensembl_no_rbh_available",
                        "orthology_type": "",
                        "percent_identity": "",
                        "coverage": "",
                        "evalue": "",
                        "HD": hd,
                        "JO": jo,
                        "note": f"ensembl_no_hits; error={ensembl_error}"
                    })
                else:
                    # RBH
                    try:
                        with tempfile.TemporaryDirectory(dir=tmpdir_global) as qtmp:
                            tgt, metrics = run_rbh_for_query(qid, aedes_index, culex_db_base, aedes_db_base, session, user_agent, qtmp)
                        if tgt:
                            out_rows.append({
                                "query_id": qid,
                                "target_gene_id": tgt,
                                "method_used": metrics.get("method_used", "rbh"),
                                "orthology_type": "",
                                "percent_identity": metrics.get("pident", ""),
                                "coverage": metrics.get("coverage", ""),
                                "evalue": metrics.get("evalue", ""),
                                "HD": hd,
                                "JO": jo,
                                "note": ""
                            })
                        else:
                            out_rows.append({
                                "query_id": qid,
                                "target_gene_id": "",
                                "method_used": metrics.get("method_used", "rbh_failed"),
                                "orthology_type": "",
                                "percent_identity": metrics.get("pident", ""),
                                "coverage": metrics.get("coverage", ""),
                                "evalue": metrics.get("evalue", ""),
                                "HD": hd,
                                "JO": jo,
                                "note": metrics.get("note", "")
                            })
                    except Exception as e:
                        out_rows.append({
                            "query_id": qid,
                            "target_gene_id": "",
                            "method_used": "rbh_exception",
                            "orthology_type": "",
                            "percent_identity": "",
                            "coverage": "",
                            "evalue": "",
                            "HD": hd,
                            "JO": jo,
                            "note": str(e)
                        })
            time.sleep(sleep)

        # write output CSV
        fieldnames = ["query_id","target_gene_id","method_used","orthology_type",
                      "percent_identity","coverage","evalue","HD","JO","note"]
        with open(output_csv, "w", newline="") as outf:
            writer = csv.DictWriter(outf, fieldnames=fieldnames)
            writer.writeheader()
            for r in out_rows:
                writer.writerow(r)
        return output_csv, len(out_rows)
    finally:
        # keep tmpdir for debugging if you want; here we remove it
        try:
            shutil.rmtree(tmpdir_global)
        except Exception:
            pass

# ==== GUI runner ====
def run_gui():
    root = tk.Tk()
    root.withdraw()
    messagebox.showinfo("Instructions", "Select input CSV (with GENE_ID, HD, JO). For RBH fallback you will be asked to supply local proteomes (recommended).")
    input_csv = filedialog.askopenfilename(title="Select input CSV", filetypes=[("CSV files","*.csv"),("All files","*.*")])
    if not input_csv:
        messagebox.showwarning("Canceled", "No input selected. Exiting.")
        return
    output_csv = filedialog.asksaveasfilename(title="Save output CSV as", defaultextension=".csv",
                                              filetypes=[("CSV files","*.csv"),("All files","*.*")],
                                              initialfile="aael_to_culex_mapping_with_rbh.csv")
    if not output_csv:
        messagebox.showwarning("Canceled", "No output selected. Exiting.")
        return
    email = None
    while not email:
        email = simpledialog.askstring("Contact email", "Enter contact email (required for polite API access):")
        if email is None:
            messagebox.showwarning("Canceled", "Email is required. Exiting.")
            return
        email = email.strip()
    target_species = simpledialog.askstring("Target species substring", "Enter substring to identify target species in Ensembl (default 'culex'):", initialvalue="culex")
    if target_species is None:
        target_species = "culex"
    # ask for local Culex proteome FASTA (required for RBH)
    use_culex = messagebox.askyesno("Culex proteome", "Do you have a local Culex proteome FASTA to use for RBH? (recommended)", default=messagebox.YES)
    culex_fasta = None
    if use_culex:
        culex_fasta = filedialog.askopenfilename(title="Select Culex proteome FASTA", filetypes=[("FASTA","*.fa *.faa *.fasta"),("All files","*.*")])
        if not culex_fasta:
            messagebox.showwarning("No file", "No Culex FASTA selected. RBH will not be available.")
            culex_fasta = None
    # ask for local Aedes proteome FASTA (optional)
    have_aedes = messagebox.askyesno("Aedes proteome", "Do you have a local Aedes proteome FASTA? (optional but recommended)", default=messagebox.NO)
    aedes_fasta = None
    if have_aedes:
        aedes_fasta = filedialog.askopenfilename(title="Select Aedes proteome FASTA", filetypes=[("FASTA","*.fa *.faa *.fasta"),("All files","*.*")])
        if not aedes_fasta:
            aedes_fasta = None
    # thresholds
    try:
        idmin = simpledialog.askfloat("BLAST pct identity", "Minimum percent identity to accept BLAST hits (default 30.0):", initialvalue=BLAST_IDENTITY_MIN)
        if idmin is not None:
            global BLAST_IDENTITY_MIN
            BLAST_IDENTITY_MIN = float(idmin)
    except Exception:
        pass
    try:
        covmin = simpledialog.askfloat("BLAST coverage", "Minimum fraction coverage of query (0-1) (default 0.5):", initialvalue=BLAST_COVERAGE_MIN)
        if covmin is not None:
            global BLAST_COVERAGE_MIN
            BLAST_COVERAGE_MIN = float(covmin)
    except Exception:
        pass
    try:
        sleep = simpledialog.askfloat("Sleep", "Seconds to sleep between web/API requests (default 0.25):", initialvalue=0.25)
        if sleep is None:
            sleep = 0.25
    except Exception:
        sleep = 0.25

    try:
        out_path, nrows = map_with_ensembl_and_rbh(input_csv, output_csv, email, target_species, culex_fasta, aedes_fasta, sleep=sleep)
        messagebox.showinfo("Done", f"Mapping complete. Wrote {nrows} rows to:\n{out_path}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred:\n{e}")

if __name__ == "__main__":
    run_gui()
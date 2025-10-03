#!/usr/bin/env python3
"""
Map Aedes aegypti AAEL gene IDs to Culex pipiens pallens orthologs using:
  - Ensembl/VectorBase homology (primary)
  - Reciprocal-BLAST (RBH) fallback (high accuracy)

This variant can automatically download proteomes from NCBI (RefSeq preferred, GenBank fallback)
and build BLAST DBs before running RBH.

Usage:
    python ensembl_map_aael_to_culex_rbh_ncbi.py

Dependencies:
    pip install pandas requests biopython
    BLAST+ (makeblastdb, blastp) installed and on PATH
"""
import os
import sys
import csv
import time
import gzip
import shutil
import subprocess
import tempfile
from io import StringIO
from typing import Optional, Dict, Tuple, List

import requests
import pandas as pd
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord

# ==== Config defaults ====
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS_JSON = {"Accept": "application/json"}
HEADERS_FASTA = {"Accept": "text/x-fasta"}
BLAST_EVALUE = 1e-5
BLAST_IDENTITY_MIN = 30.0
BLAST_COVERAGE_MIN = 0.5
BLAST_MAX_TARGETS = 5

# ==== Helper utilities ====
def which(exe_name: str) -> Optional[str]:
    return shutil.which(exe_name)

def check_blast_tools():
    for exe in ("makeblastdb", "blastp"):
        if which(exe) is None:
            raise FileNotFoundError(f"{exe} not found on PATH. Install BLAST+ and ensure {exe} is available.")

def make_blast_db(fasta_path: str, db_base: str):
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
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"blastp failed: {proc.stderr.strip()}")
    return proc.stdout

# ==== Ensembl helpers ====
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

# ==== NCBI proteome download helpers ====
def ncbi_find_assembly_ftp(species_name: str, email: str, prefer_refseq: bool = True) -> Optional[Dict]:
    """
    Find an assembly summary via Entrez/assembly and return a dict containing ftp paths.
    Returns {'refseq_ftp':..., 'genbank_ftp':..., 'assembly_accession':...} or None
    """
    Entrez.email = email
    # search assembly db for organism
    term = f'{species_name}[Organism]'
    try:
        handle = Entrez.esearch(db="assembly", term=term, retmax=20)
        res = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise RuntimeError(f"Entrez esearch failed: {e}")

    idlist = res.get("IdList", [])
    if not idlist:
        return None

    # prefer assemblies with 'refseq' in summary if requested; otherwise take first
    # fetch summaries
    try:
        sumh = Entrez.esummary(db="assembly", id=",".join(idlist))
        summ = Entrez.read(sumh)
        sumh.close()
    except Exception as e:
        raise RuntimeError(f"Entrez esummary failed: {e}")

    # structure: summ['DocumentSummarySet']['DocumentSummary'] list; but Entrez.read may flatten
    doclist = summ.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    if not doclist:
        # fallback: use raw summ as list/dict
        if isinstance(summ, list):
            doclist = summ
        else:
            doclist = []

    # choose best assembly: prefer RefSeq, then latest (highest version)
    chosen = None
    for doc in doclist:
        refftp = doc.get("FtpPath_RefSeq") or ""
        genftp = doc.get("FtpPath_GenBank") or ""
        asm_accession = doc.get("AssemblyAccession") or doc.get("assembly_accession") or ""
        if prefer_refseq and refftp:
            chosen = {"refseq_ftp": refftp, "genbank_ftp": genftp, "assembly_accession": asm_accession}
            break
        if not chosen and (refftp or genftp):
            chosen = {"refseq_ftp": refftp, "genbank_ftp": genftp, "assembly_accession": asm_accession}
    return chosen

def try_download_protein_fasta_from_ftp(ftp_path: str, dest_dir: str) -> Optional[str]:
    """
    Given ftp_path like 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../GCF_XXX.1', try to download typical
    protein FASTA filenames and return path to local fasta (gunzipped if needed), or None on failure.
    """
    if ftp_path.startswith("ftp://"):
        http_base = ftp_path.replace("ftp://", "https://")
    elif ftp_path.startswith("https://"):
        http_base = ftp_path
    else:
        http_base = "https://" + ftp_path

    # basename is last segment of ftp_path
    base = os.path.basename(ftp_path.rstrip("/"))
    # candidate filenames
    candidates = [
        f"{base}_protein.faa.gz",
        f"{base}_protein.faa",
        f"{base}_protein.faa.gz",  # repeated for clarity
        "protein.faa.gz",
        "protein.faa",
    ]
    session = requests.Session()
    for c in candidates:
        url = http_base + "/" + c
        try:
            # HEAD first
            h = session.head(url, timeout=20)
            if h.status_code == 200:
                # download
                local = os.path.join(dest_dir, c)
                with session.get(url, stream=True, timeout=60) as r:
                    r.raise_for_status()
                    with open(local, "wb") as fh:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                fh.write(chunk)
                # if gz, gunzip
                if local.endswith(".gz"):
                    out_path = local[:-3]
                    with gzip.open(local, "rb") as fin, open(out_path, "wb") as fout:
                        shutil.copyfileobj(fin, fout)
                    os.remove(local)
                    return out_path
                else:
                    return local
        except Exception:
            continue
    return None

# ==== FASTA helpers ====
def index_fasta_by_id_substring(fasta_path: str) -> Dict[str, SeqRecord]:
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

# ==== RBH logic (uses local DBs) ====
def run_rbh_for_query(query_gene: str, aedes_index: Optional[Dict[str, SeqRecord]],
                      culex_db_base: str, aedes_db_base: Optional[str],
                      session: requests.Session, user_agent: str,
                      tmpdir: str) -> Tuple[Optional[str], Dict]:
    # obtain query protein sequence
    seqrec = None
    if aedes_index is not None:
        matches = find_records_containing_substr(aedes_index, query_gene)
        if matches:
            seqrec = matches[0]
    if seqrec is None:
        seqrec = ensembl_fetch_protein_sequence(query_gene, session, user_agent)
    if seqrec is None:
        return None, {"method_used": "rbh_no_query_seq", "note": "no_query_sequence"}
    qf = os.path.join(tmpdir, f"query_{query_gene}.fa")
    SeqIO.write(seqrec, qf, "fasta")
    out = run_blastp_query(qf, culex_db_base)
    if not out.strip():
        return None, {"method_used": "rbh_no_hit", "note": "no_hit_to_culex"}
    top_line = out.strip().splitlines()[0]
    parts = top_line.split('\t')
    sseqid = parts[1]
    pident = float(parts[2])
    align_len = int(parts[3])
    evalue = float(parts[10])
    qlen = len(seqrec.seq)
    coverage = align_len / qlen if qlen > 0 else 0.0
    if pident < BLAST_IDENTITY_MIN or coverage < BLAST_COVERAGE_MIN or evalue > BLAST_EVALUE:
        return None, {"method_used": "rbh_hit_low_quality", "pident": pident, "coverage": coverage, "evalue": evalue, "note": "hit_below_thresholds"}
    # retrieve subject seq using blastdbcmd if available, else try to read from recorded source fasta
    subj_fa = os.path.join(tmpdir, f"subj_{sseqid}.fa")
    if which("blastdbcmd"):
        try:
            cmd = ["blastdbcmd", "-db", culex_db_base, "-entry", sseqid, "-outfmt", "%f"]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            SeqIO.write(list(SeqIO.parse(StringIO(res.stdout), "fasta")), subj_fa, "fasta")
        except Exception:
            return None, {"method_used": "rbh_no_subject_seq", "note": "cannot_retrieve_subject_seq"}
    else:
        return None, {"method_used": "rbh_no_subject_seq_no_blastdbcmd", "note": "blastdbcmd_missing"}
    # blast subject -> aedes DB
    aedes_db_for_rbh = aedes_db_base
    if aedes_db_for_rbh is None:
        if aedes_index is None:
            return None, {"method_used": "rbh_no_aedes_db", "note": "no_aedes_db_available"}
        aedes_tmp_fa = os.path.join(tmpdir, "aedes_from_index.fa")
        SeqIO.write(list(aedes_index.values()), aedes_tmp_fa, "fasta")
        aedes_db_tmp = os.path.join(tmpdir, "aedes_tmp_db")
        make_blast_db(aedes_tmp_fa, aedes_db_tmp)
        aedes_db_for_rbh = aedes_db_tmp
    out2 = run_blastp_query(subj_fa, aedes_db_for_rbh)
    if not out2.strip():
        return None, {"method_used": "rbh_no_reciprocal_hit", "note": "no_reciprocal_hit"}
    top2 = out2.strip().splitlines()[0].split('\t')
    reciprocal_hit_id = top2[1]
    # accept reciprocal if query_gene appears in reciprocal_hit_id or qid
    if query_gene in reciprocal_hit_id or query_gene in top2[0]:
        return reciprocal_hit_id, {"method_used": "rbh", "pident": pident, "coverage": coverage, "evalue": evalue}
    else:
        return None, {"method_used": "rbh_nonreciprocal", "pident": pident, "coverage": coverage, "evalue": evalue, "reciprocal_top": reciprocal_hit_id}

# ==== Main flow combining Ensembl + RBH + NCBI download ====
def map_with_ensembl_and_rbh_ncbi(input_csv: str, output_csv: str, email: str,
                                  target_species: str, auto_download_ncbi: bool,
                                  sleep: float = 0.25):
    check_blast_tools()
    session = requests.Session()
    user_agent = f"ensembl-rbh-ncbi ({email})" if email else "ensembl-rbh-ncbi"
    df = pd.read_csv(input_csv)
    if "GENE_ID" not in df.columns:
        raise ValueError("Input CSV must contain GENE_ID column")
    tmpdir_global = tempfile.mkdtemp(prefix="ensembl_rbh_ncbi_")
    culex_db_base = None
    aedes_db_base = None
    aedes_index = None

    try:
        # If auto download requested, attempt to download Culex and Aedes proteomes from NCBI
        culex_fasta = None
        aedes_fasta = None
        if auto_download_ncbi:
            # prompt species names (use target_species for Culex)
            # target_species param likely 'culex' substring; but for NCBI we need full organism name
            culex_species_name = simpledialog.askstring("NCBI species name (culex)", "Enter NCBI organism name for Culex (example: 'Culex pipiens pallens'):", initialvalue="Culex pipiens pallens")
            if culex_species_name:
                try:
                    asm = ncbi_find_assembly_ftp(culex_species_name, email, prefer_refseq=True)
                    if asm is not None:
                        ftp = asm.get("refseq_ftp") or asm.get("genbank_ftp")
                        if ftp:
                            messagebox.showinfo("Downloading", f"Attempting to download Culex proteome from:\n{ftp}")
                            culex_fasta = try_download_protein_fasta_from_ftp(ftp, tmpdir_global)
                            if culex_fasta:
                                # record source path
                                culex_db_base = os.path.join(tmpdir_global, "culex_db")
                                make_blast_db(culex_fasta, culex_db_base)
                                # store source for later retrieval attempts
                                with open(culex_db_base + ".source_fasta", "w") as fh:
                                    fh.write(os.path.abspath(culex_fasta))
                            else:
                                messagebox.showwarning("Not found", f"Could not find protein FASTA at {ftp}. You will be asked to provide a local file.")
                    else:
                        messagebox.showwarning("No assembly", f"No assemblies found for '{culex_species_name}' on NCBI.")
                except Exception as e:
                    messagebox.showwarning("NCBI error", f"Error while querying NCBI for Culex: {e}")
            # Aedes
            aedes_species_name = simpledialog.askstring("NCBI species name (aedes)", "Enter NCBI organism name for Aedes (example: 'Aedes aegypti'):", initialvalue="Aedes aegypti")
            if aedes_species_name:
                try:
                    asm2 = ncbi_find_assembly_ftp(aedes_species_name, email, prefer_refseq=True)
                    if asm2 is not None:
                        ftp2 = asm2.get("refseq_ftp") or asm2.get("genbank_ftp")
                        if ftp2:
                            messagebox.showinfo("Downloading", f"Attempting to download Aedes proteome from:\n{ftp2}")
                            aedes_fasta = try_download_protein_fasta_from_ftp(ftp2, tmpdir_global)
                            if aedes_fasta:
                                aedes_db_base = os.path.join(tmpdir_global, "aedes_db")
                                make_blast_db(aedes_fasta, aedes_db_base)
                                aedes_index = index_fasta_by_id_substring(aedes_fasta)
                            else:
                                messagebox.showwarning("Not found", f"Could not find protein FASTA at {ftp2}. You will be asked to provide a local file.")
                    else:
                        messagebox.showwarning("No assembly", f"No assemblies found for '{aedes_species_name}' on NCBI.")
                except Exception as e:
                    messagebox.showwarning("NCBI error", f"Error while querying NCBI for Aedes: {e}")

        # If any DB missing, ask user to provide local FASTAs
        if culex_db_base is None:
            use_local = messagebox.askyesno("Local Culex FASTA", "Provide local Culex proteome FASTA for RBH? (Yes = choose file)", default=messagebox.YES)
            if use_local:
                cpath = filedialog.askopenfilename(title="Select Culex proteome FASTA", filetypes=[("FASTA","*.fa *.faa *.fasta"),("All files","*.*")])
                if cpath:
                    culex_fasta = cpath
                    culex_db_base = os.path.join(tmpdir_global, "culex_db")
                    make_blast_db(culex_fasta, culex_db_base)
                    with open(culex_db_base + ".source_fasta", "w") as fh:
                        fh.write(os.path.abspath(culex_fasta))
        if aedes_db_base is None:
            use_local2 = messagebox.askyesno("Local Aedes FASTA", "Provide local Aedes proteome FASTA? (recommended; Yes = choose file)", default=messagebox.NO)
            if use_local2:
                apath = filedialog.askopenfilename(title="Select Aedes proteome FASTA", filetypes=[("FASTA","*.fa *.faa *.fasta"),("All files","*.*")])
                if apath:
                    aedes_fasta = apath
                    aedes_db_base = os.path.join(tmpdir_global, "aedes_db")
                    make_blast_db(aedes_fasta, aedes_db_base)
                    aedes_index = index_fasta_by_id_substring(aedes_fasta)

        # Now iterate input and perform Ensembl + RBH mapping
        out_rows = []
        total = len(df)
        for i, row in df.iterrows():
            qid = str(row["GENE_ID"]).strip().strip('"')
            hd = row.get("HD", "")
            jo = row.get("JO", "")
            # Ensembl attempt
            try:
                homs = ensembl_homology(qid, requests.Session(), user_agent, target_species)
            except Exception as e:
                homs = []
                ensembl_error = str(e)
            else:
                ensembl_error = ""
            chosen = []
            if homs:
                for h in homs:
                    typ = h.get("orthology_type","").lower()
                    if "one2one" in typ or "one_to_one" in typ or "one-to-one" in typ:
                        chosen.append((h,"ensembl_one2one"))
                if not chosen:
                    for h in homs:
                        chosen.append((h,"ensembl_other"))
            if chosen:
                for h,tag in chosen:
                    out_rows.append({
                        "query_id": qid,
                        "target_gene_id": h.get("target_gene_id",""),
                        "method_used": "ensembl",
                        "orthology_type": h.get("orthology_type",""),
                        "percent_identity": h.get("perc_id",""),
                        "coverage": "",
                        "evalue": "",
                        "HD": hd,
                        "JO": jo,
                        "note": tag
                    })
            else:
                # RBH fallback
                if culex_db_base is None:
                    out_rows.append({
                        "query_id": qid,
                        "target_gene_id": "",
                        "method_used": "no_rbh_available",
                        "orthology_type": "",
                        "percent_identity": "",
                        "coverage": "",
                        "evalue": "",
                        "HD": hd,
                        "JO": jo,
                        "note": f"ensembl_no_hits; error={ensembl_error}"
                    })
                else:
                    try:
                        with tempfile.TemporaryDirectory(dir=tmpdir_global) as qtmp:
                            tgt, metrics = run_rbh_for_query(qid, aedes_index, culex_db_base, aedes_db_base, requests.Session(), user_agent, qtmp)
                        if tgt:
                            out_rows.append({
                                "query_id": qid,
                                "target_gene_id": tgt,
                                "method_used": metrics.get("method_used","rbh"),
                                "orthology_type": "",
                                "percent_identity": metrics.get("pident",""),
                                "coverage": metrics.get("coverage",""),
                                "evalue": metrics.get("evalue",""),
                                "HD": hd,
                                "JO": jo,
                                "note": ""
                            })
                        else:
                            out_rows.append({
                                "query_id": qid,
                                "target_gene_id": "",
                                "method_used": metrics.get("method_used","rbh_failed"),
                                "orthology_type": "",
                                "percent_identity": metrics.get("pident",""),
                                "coverage": metrics.get("coverage",""),
                                "evalue": metrics.get("evalue",""),
                                "HD": hd,
                                "JO": jo,
                                "note": metrics.get("note","")
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

        # write CSV
        fieldnames = ["query_id","target_gene_id","method_used","orthology_type",
                      "percent_identity","coverage","evalue","HD","JO","note"]
        with open(output_csv, "w", newline="") as outf:
            writer = csv.DictWriter(outf, fieldnames=fieldnames)
            writer.writeheader()
            for r in out_rows:
                writer.writerow(r)
        return output_csv, len(out_rows)
    finally:
        # remove tempdir for cleanliness
        try:
            shutil.rmtree(tmpdir_global)
        except Exception:
            pass

# ==== GUI ====
def run_gui():
    root = tk.Tk()
    root.withdraw()
    messagebox.showinfo("Instructions", "Select input CSV (with GENE_ID, HD, JO). You can choose automatic NCBI download for proteomes.")
    input_csv = filedialog.askopenfilename(title="Select input CSV", filetypes=[("CSV files","*.csv"),("All files","*.*")])
    if not input_csv:
        messagebox.showwarning("Canceled", "No input selected. Exiting.")
        return
    output_csv = filedialog.asksaveasfilename(title="Save output CSV as", defaultextension=".csv",
                                              filetypes=[("CSV files","*.csv"),("All files","*.*")],
                                              initialfile="aael_to_culex_mapping_with_rbh_ncbi.csv")
    if not output_csv:
        messagebox.showwarning("Canceled", "No output selected. Exiting.")
        return
    email = None
    while not email:
        email = simpledialog.askstring("Contact email", "Enter contact email (required for Entrez and polite API use):")
        if email is None:
            messagebox.showwarning("Canceled", "Email is required. Exiting.")
            return
        email = email.strip()
    target_species = simpledialog.askstring("Target species substring", "Enter substring to identify target species in Ensembl (default 'culex'):", initialvalue="culex")
    if target_species is None:
        target_species = "culex"
    auto = messagebox.askyesno("Automatic NCBI download", "Attempt to download proteomes automatically from NCBI (RefSeq preferred)?", default=messagebox.YES)
    try:
        out_path, nrows = map_with_ensembl_and_rbh_ncbi(input_csv, output_csv, email, target_species, auto)
        messagebox.showinfo("Done", f"Mapping complete. Wrote {nrows} rows to:\n{out_path}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred:\n{e}")

if __name__ == "__main__":
    run_gui()
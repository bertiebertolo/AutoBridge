#!/usr/bin/env python3
"""
BLAST+ helper functions used by AUTOBridge GUI and RBH workflow.
"""
import os
import subprocess
import shutil

BLAST_MAX_TARGETS = 5
BLAST_EVALUE = 1e-5

def which(exe_name: str):
    return shutil.which(exe_name)

def check_blast_tools():
    for exe in ("makeblastdb", "blastp", "blastn", "blastdbcmd"):
        if which(exe) is None:
            raise FileNotFoundError(f"{exe} not found on PATH. Install BLAST+ and ensure {exe} is available.")

def make_blast_db(fasta_path: str, db_base: str, dbtype: str = "prot"):
    if os.path.exists(db_base + ".pin") or os.path.exists(db_base + ".psq") or os.path.exists(db_base + ".phr"):
        return db_base
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", dbtype, "-out", db_base]
    subprocess.run(cmd, check=True)
    return db_base

def run_blastp(query_fasta: str, db_base: str, out_path: str, outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
    cmd = [
        "blastp", "-query", query_fasta, "-db", db_base,
        "-outfmt", outfmt, "-max_target_seqs", str(BLAST_MAX_TARGETS),
        "-evalue", str(BLAST_EVALUE), "-out", out_path
    ]
    subprocess.run(cmd, check=True)
    return out_path

def run_blastn(query_fasta: str, db_base: str, out_path: str, outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
    cmd = [
        "blastn", "-query", query_fasta, "-db", db_base,
        "-outfmt", outfmt, "-max_target_seqs", str(BLAST_MAX_TARGETS),
        "-evalue", str(BLAST_EVALUE), "-out", out_path
    ]
    subprocess.run(cmd, check=True)
    return out_path

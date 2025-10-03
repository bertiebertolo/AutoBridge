#!/usr/bin/env python3
"""
Helpers to fetch proteomes from NCBI.

Primary entrypoint:
  get_proteome_for_taxon(taxon: str, email: str, out_dir: Optional[str] = None) -> Optional[str]

Behavior:
- First tries `ncbi-datasets-cli` (recommended) if available on PATH.
- If datasets is not available, falls back to Entrez/Biopython to locate an assembly FTP
  and attempts to download a *_protein.faa(.gz) file.
- Returns path to an uncompressed protein FASTA file, or None if no proteome could be obtained.
"""
from typing import Optional
import os
import shutil
import subprocess
import glob
import gzip
import requests

# Biopython Entrez optional
try:
    from Bio import Entrez
    BIOPYTHON_AVAILABLE = True
except Exception:
    BIOPYTHON_AVAILABLE = False


def _run_cmd(cmd: str, cwd: Optional[str] = None) -> str:
    """Run shell command, raise on non-zero exit, return stdout."""
    p = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\nstdout: {p.stdout}\nstderr: {p.stderr}")
    return p.stdout


def _find_protein_fasta_in_extracted(path: str) -> Optional[str]:
    """Search extracted dataset dir for likely protein FASTA files."""
    patterns = [
        "**/*protein*.faa*",
        "**/*protein*.fasta*",
        "**/*proteins*.faa*",
        "**/*proteins*.fasta*",
    ]
    for pat in patterns:
        matches = glob.glob(os.path.join(path, pat), recursive=True)
        if matches:
            return matches[0]
    # fallback: any .faa/.fasta
    for pat in ["**/*.faa*", "**/*.fasta*"]:
        matches = glob.glob(os.path.join(path, pat), recursive=True)
        if matches:
            return matches[0]
    return None


def _uncompress_if_gz(path: str) -> str:
    """If path ends with .gz, gunzip to same dir and return new path; otherwise return original."""
    if not path.endswith(".gz"):
        return path
    outpath = os.path.splitext(path)[0]
    with gzip.open(path, "rb") as f_in, open(outpath, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return outpath


def download_proteome_via_datasets(taxon: str, out_dir: Optional[str] = None) -> Optional[str]:
    """
    Use `datasets` CLI to download protein FASTA for a taxon.
    Returns path to an uncompressed FASTA or None.
    """
    if shutil.which("datasets") is None:
        return None
    if out_dir is None:
        out_dir = os.getcwd()
    # safe filename and extraction folder inside out_dir
    safe = taxon.replace(" ", "_").replace("/", "_")
    zipname = os.path.join(out_dir, f"datasets_{safe}.zip")
    try:
        cmd = f'datasets download genome taxon "{taxon}" --filename "{zipname}" --include protein'
        _run_cmd(cmd, cwd=out_dir)
        extract_dir = os.path.join(out_dir, f"datasets_{safe}_extracted")
        os.makedirs(extract_dir, exist_ok=True)
        shutil.unpack_archive(zipname, extract_dir)
        pf = _find_protein_fasta_in_extracted(extract_dir)
        if pf is None:
            return None
        return _uncompress_if_gz(pf)
    except Exception:
        # cleanup partial files if any
        try:
            if os.path.exists(zipname):
                os.remove(zipname)
        except Exception:
            pass
        return None


def download_proteome_via_entrez(taxon: str, email: str, out_dir: Optional[str] = None) -> Optional[str]:
    """
    Entrez fallback: search assembly db for organism, get FTP path(s), and try to download proteins.
    Requires Biopython Entrez (BIOPYTHON_AVAILABLE).
    Returns path to an uncompressed FASTA or None.
    """
    if not BIOPYTHON_AVAILABLE:
        return None
    if out_dir is None:
        out_dir = os.getcwd()
    Entrez.email = email
    try:
        handle = Entrez.esearch(db="assembly", term=f"{taxon}[Organism]", retmax=25)
        rec = Entrez.read(handle)
        handle.close()
        ids = rec.get("IdList", [])
        if not ids:
            return None
    except Exception:
        return None

    for aid in ids:
        try:
            h = Entrez.esummary(db="assembly", id=aid, report="full")
            summary = Entrez.read(h)
            h.close()
            docs = summary.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            for d in docs:
                ftp = d.get("FtpPath_RefSeq") or d.get("FtpPath_GenBank")
                if not ftp:
                    continue
                base = os.path.basename(ftp)
                # common filenames to try
                candidates = [
                    f"{ftp}/{base}_protein.faa.gz",
                    f"{ftp}/{base}_protein.faa",
                    f"{ftp}/{base}_translated_cds.faa.gz",
                    f"{ftp}/{base}_protein.faa.gz"
                ]
                # try HEAD for candidates
                found_url = None
                for c in candidates:
                    try:
                        rr = requests.head(c, timeout=10)
                        if rr.status_code == 200:
                            found_url = c
                            break
                    except Exception:
                        continue
                # if not found by HEAD, try to scan directory HTML (best-effort)
                if not found_url:
                    try:
                        r = requests.get(ftp + "/", timeout=20)
                        if r.status_code == 200:
                            import re
                            m = re.search(r'href="([^"]*protein[^"]*\.faa(\.gz)?)"', r.text, flags=re.IGNORECASE)
                            if m:
                                found_url = ftp + "/" + m.group(1)
                    except Exception:
                        pass
                if not found_url:
                    continue
                # download
                local_fname = os.path.join(out_dir, os.path.basename(found_url))
                with requests.get(found_url, stream=True, timeout=60) as rr:
                    if rr.status_code != 200:
                        continue
                    with open(local_fname, "wb") as fh:
                        for chunk in rr.iter_content(chunk_size=8192):
                            if chunk:
                                fh.write(chunk)
                return _uncompress_if_gz(local_fname)
        except Exception:
            continue
    return None


def get_proteome_for_taxon(taxon: str, email: str, out_dir: Optional[str] = None) -> Optional[str]:
    """
    High-level helper. Tries datasets CLI first, then Entrez fallback.
    Returns path to an uncompressed protein FASTA, or None if not available.
    """
    # Try datasets CLI
    try:
        pf = download_proteome_via_datasets(taxon, out_dir=out_dir)
        if pf:
            return pf
    except Exception:
        pass
    # Entrez fallback
    try:
        pf = download_proteome_via_entrez(taxon, email, out_dir=out_dir)
        if pf:
            return pf
    except Exception:
        pass
    return None

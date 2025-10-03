#!/usr/bin/env python3
"""
Minimal NCBI assembly/proteome helpers using Biopython Entrez; used by AUTOBridge if requested.

Functions:
 - ncbi_find_assembly_ftp(species_name, email)
 - try_download_protein_fasta_from_ftp(ftp_path, dest_dir)
"""
import os
import gzip
import shutil
import requests
from Bio import Entrez
from typing import Optional

def ncbi_find_assembly_ftp(species_name: str, email: str, prefer_refseq: bool = True):
    Entrez.email = email
    handle = Entrez.esearch(db="assembly", term=f"{species_name}[Organism]", retmax=20)
    res = Entrez.read(handle)
    handle.close()
    ids = res.get("IdList", [])
    if not ids:
        return None
    sumh = Entrez.esummary(db="assembly", id=",".join(ids))
    summ = Entrez.read(sumh)
    sumh.close()
    doclist = summ.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    chosen = None
    for doc in doclist:
        refftp = doc.get("FtpPath_RefSeq") or ""
        genftp = doc.get("FtpPath_GenBank") or ""
        asm_accession = doc.get("AssemblyAccession") or ""
        if prefer_refseq and refftp:
            chosen = {"refseq_ftp": refftp, "genbank_ftp": genftp, "assembly_accession": asm_accession}
            break
        if not chosen and (refftp or genftp):
            chosen = {"refseq_ftp": refftp, "genbank_ftp": genftp, "assembly_accession": asm_accession}
    return chosen

def try_download_protein_fasta_from_ftp(ftp_path: str, dest_dir: str):
    if ftp_path.startswith("ftp://"):
        http_base = ftp_path.replace("ftp://", "https://")
    elif ftp_path.startswith("https://"):
        http_base = ftp_path
    else:
        http_base = "https://" + ftp_path
    base = os.path.basename(ftp_path.rstrip("/"))
    candidates = [f"{base}_protein.faa.gz", f"{base}_protein.faa", "protein.faa.gz", "protein.faa"]
    session = requests.Session()
    for c in candidates:
        url = http_base + "/" + c
        try:
            h = session.head(url, timeout=20)
            if h.status_code == 200:
                local = os.path.join(dest_dir, c)
                with session.get(url, stream=True, timeout=60) as r:
                    r.raise_for_status()
                    with open(local, "wb") as fh:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                fh.write(chunk)
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
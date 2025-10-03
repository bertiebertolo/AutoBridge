#!/usr/bin/env python3
"""
Ensembl homology + RBH orchestration functions used by AUTOBridge.

Core functionality:
 - Ensembl homology lookup (REST)
 - RBH fallback using BLAST+ and local proteome FASTAs (or Ensembl sequences)
"""
import os
import tempfile
import requests
import time
import csv
from typing import Optional
from Bio import SeqIO
from io import StringIO

from blast_utils import make_blast_db, run_blastp, check_blast_tools

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS_JSON = {"Accept": "application/json"}
HEADERS_FASTA = {"Accept": "text/x-fasta"}

def ensembl_homology(gene_id: str, target_species_substr: str, user_agent: str):
    url = f"{ENSEMBL_REST}/homology/id/{gene_id}"
    params = {"type": "orthologues", "sequence": 0, "format": "full"}
    headers = HEADERS_JSON.copy()
    headers["User-Agent"] = user_agent
    r = requests.get(url, params=params, headers=headers, timeout=30)
    r.raise_for_status()
    data = r.json()
    out = []
    for rec in data.get("data", []):
        for h in rec.get("homologies", []):
            tgt = h.get("target", {})
            species = tgt.get("species","")
            if target_species_substr.lower() not in species.lower():
                continue
            out.append({
                "target_gene_id": tgt.get("gene_id") or tgt.get("id",""),
                "target_protein_id": tgt.get("id",""),
                "target_species": species,
                "orthology_type": h.get("type",""),
                "confidence": h.get("confidence",""),
                "perc_id": h.get("percentage_identities") or h.get("perc_id","")
            })
    return out

def ensembl_fetch_protein_sequence(identifier: str, user_agent: str):
    url = f"{ENSEMBL_REST}/sequence/id/{identifier}"
    params = {"type": "protein"}
    headers = HEADERS_FASTA.copy()
    headers["User-Agent"] = user_agent
    r = requests.get(url, params=params, headers=headers, timeout=30)
    if r.status_code == 200:
        fasta = r.text
        try:
            rec = next(SeqIO.parse(StringIO(fasta), "fasta"))
            return rec
        except Exception:
            return None
    return None

def map_genes(input_csv: str, output_csv: str, email: str, target_species: str, culex_db_fasta: Optional[str]=None, aedes_db_fasta: Optional[str]=None, sleep: float=0.25):
    """
    Batch mapping entrypoint used by AUTOBridge (non-NCBI-auto variant).
    """
    check_blast_tools()
    df_in = []
    with open(input_csv) as inf:
        reader = csv.DictReader(inf)
        for r in reader:
            df_in.append(r)
    aedes_index = None
    if aedes_db_fasta:
        aedes_index = {rec.id:rec for rec in SeqIO.parse(aedes_db_fasta, "fasta")}
        aedes_db_base = os.path.splitext(aedes_db_fasta)[0] + "_db"
        make_blast_db(aedes_db_fasta, aedes_db_base)
    else:
        aedes_db_base = None
    if culex_db_fasta:
        culex_db_base = os.path.splitext(culex_db_fasta)[0] + "_db"
        make_blast_db(culex_db_fasta, culex_db_base)
    else:
        culex_db_base = None

    user_agent = f"AUTOBridge ({email})"
    out_rows = []
    for row in df_in:
        qid = row.get("GENE_ID")
        hd = row.get("HD","")
        jo = row.get("JO","")
        try:
            homs = ensembl_homology(qid, target_species, user_agent)
        except Exception as e:
            homs = []
            ensembl_err = str(e)
        else:
            ensembl_err = ""
        chosen = []
        if homs:
            for h in homs:
                typ = h.get("orthology_type","").lower()
                if "one2one" in typ or "one_to_one" in typ or "one-to-one" in typ:
                    chosen.append(h)
            if not chosen:
                chosen = homs
        if chosen:
            for h in chosen:
                out_rows.append({
                    "query_id": qid,
                    "target_gene_id": h.get("target_gene_id"),
                    "method_used": "ensembl",
                    "orthology_type": h.get("orthology_type"),
                    "percent_identity": h.get("perc_id"),
                    "coverage": "",
                    "evalue": "",
                    "HD": hd,
                    "JO": jo,
                    "note": ""
                })
        else:
            # RBH fallback (simple)
            if not culex_db_base:
                out_rows.append({
                    "query_id": qid,
                    "target_gene_id": "",
                    "method_used": "ensembl_no_rbh",
                    "orthology_type": "",
                    "percent_identity": "",
                    "coverage": "",
                    "evalue": "",
                    "HD": hd,
                    "JO": jo,
                    "note": f"ensembl_no_hits; error={ensembl_err}"
                })
            else:
                seqrec = None
                if aedes_index:
                    for k,rec in aedes_index.items():
                        if qid in k or qid in rec.description:
                            seqrec = rec
                            break
                if seqrec is None:
                    seqrec = ensembl_fetch_protein_sequence(qid, user_agent)
                if seqrec is None:
                    out_rows.append({
                        "query_id": qid,
                        "target_gene_id": "",
                        "method_used": "rbh_no_query_seq",
                        "orthology_type": "",
                        "percent_identity": "",
                        "coverage": "",
                        "evalue": "",
                        "HD": hd,
                        "JO": jo,
                        "note": "no_query_sequence"
                    })
                else:
                    with tempfile.TemporaryDirectory() as td:
                        qf = os.path.join(td, f"{qid}.fa")
                        SeqIO.write(seqrec, qf, "fasta")
                        blast_out = os.path.join(td, "blastp.tsv")
                        run_blastp(qf, culex_db_base, blast_out)
                        if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
                            first = open(blast_out).read().splitlines()[0].split('\t')
                            sseqid = first[1]; pident=float(first[2]); alen=int(first[3]); evalue=float(first[10])
                            cov = alen/len(seqrec.seq) if len(seqrec.seq)>0 else 0.0
                            if pident >= 30.0 and cov >= 0.5:
                                out_rows.append({
                                    "query_id": qid,
                                    "target_gene_id": sseqid,
                                    "method_used": "rbh",
                                    "orthology_type": "",
                                    "percent_identity": pident,
                                    "coverage": cov,
                                    "evalue": evalue,
                                    "HD": hd,
                                    "JO": jo,
                                    "note": ""
                                })
                            else:
                                out_rows.append({
                                    "query_id": qid,
                                    "target_gene_id": "",
                                    "method_used": "rbh_low_quality",
                                    "orthology_type": "",
                                    "percent_identity": pident,
                                    "coverage": cov,
                                    "evalue": evalue,
                                    "HD": hd,
                                    "JO": jo,
                                    "note": "hit_below_thresholds"
                                })
                        else:
                            out_rows.append({
                                "query_id": qid,
                                "target_gene_id": "",
                                "method_used": "rbh_no_hit",
                                "orthology_type": "",
                                "percent_identity": "",
                                "coverage": "",
                                "evalue": "",
                                "HD": hd,
                                "JO": jo,
                                "note": "no_hit_to_culex"
                            })
        time.sleep(sleep)

    with open(output_csv, "w", newline="") as outf:
        fieldnames = ["query_id","target_gene_id","method_used","orthology_type","percent_identity","coverage","evalue","HD","JO","note"]
        import csv
        w = csv.DictWriter(outf, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    return output_csv, len(out_rows)
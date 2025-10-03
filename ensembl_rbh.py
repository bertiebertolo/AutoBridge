#!/usr/bin/env python3
"""
NCBI-first orthology mapping for AUTOBridge (replaces Ensembl homology flow).

This module implements a simple NCBI-datasets-driven Reciprocal Best Hit (RBH)
flow that:

- Downloads proteomes for a source taxon (e.g. "Aedes aegypti") and a target taxon
  (e.g. "Culex pipiens pallens") using ncbi_utils.get_proteome_for_taxon (datasets CLI
  preferred; Entrez fallback).
- Finds sequences in the source proteome that match input gene identifiers by
  substring matching against FASTA headers/descriptions.
- For each matched source sequence:
  - blastp the sequence against the target DB to get the top hit.
  - retrieve that target sequence and blastp it back against the source DB.
  - if the reciprocal top hit is the original source sequence, report an RBH mapping.
- Writes CSV results with mapping rows.

Notes / assumptions:
- Input CSV should have a column named GENE_ID containing identifiers (Ensembl IDs,
  gene symbols, or other identifiers that appear in the source FASTA header/description).
- This is NCBI-only: no Ensembl homology calls are made in the NCBI-RBH flow.
- Requires BLAST+ on PATH, and ncbi_utils.get_proteome_for_taxon available in repo.
"""
import os
import csv
import tempfile
import time
import logging
from typing import Optional, Tuple, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from blast_utils import check_blast_tools, make_blast_db, run_blastp, run_blastn

# ncbi_utils helper (must exist in repo, uses datasets CLI or Entrez fallback)
try:
    from ncbi_utils import get_proteome_for_taxon  # type: ignore
    NCBI_UTILS_AVAILABLE = True
except Exception:
    get_proteome_for_taxon = None  # type: ignore
    NCBI_UTILS_AVAILABLE = False

# Logging
LOGNAME = "autobridge"
DEFAULT_LOG = "autobridge.log"


def _prepare_logger(log_file: str = DEFAULT_LOG):
    logger = logging.getLogger(LOGNAME)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        fh = logging.FileHandler(log_file, mode="a")
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
        logger.addHandler(fh)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        logger.addHandler(ch)
    return logger


def _index_fasta_by_id_and_desc(fasta_path: str) -> Dict[str, SeqRecord]:
    """Return dict mapping record.id -> SeqRecord for FASTA file."""
    return {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}


def _find_source_records_for_gene_ids(source_index: Dict[str, SeqRecord], gene_ids):
    """
    For each gene_id in the input list, try to find a matching SeqRecord in source_index.
    Matching strategy: exact id match, substring in id, substring in description.
    Returns dict gene_id -> SeqRecord (or None if not found).
    """
    out = {}
    for gid in gene_ids:
        if gid is None:
            out[gid] = None
            continue
        found = None
        # exact id key
        if gid in source_index:
            found = source_index[gid]
        else:
            # substring matches in id or description
            for rid, rec in source_index.items():
                if gid in rid or gid in (rec.description or ""):
                    found = rec
                    break
        out[gid] = found
    return out


def _blast_top_hit(query_fasta_path: str, db_base: str, out_tsv: str, evalue_thresh: float = 1e-5):
    """
    Run blastp of query_fasta_path against db_base, write tabular output to out_tsv.
    Returns parsed top hit as tuple (sseqid, pident, alen, evalue) or None if no hits.
    """
    run_blastp(query_fasta_path, db_base, out_tsv)
    if not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0:
        return None
    first = open(out_tsv).read().splitlines()[0].split('\t')
    sseqid = first[1]
    pident = float(first[2])
    alen = int(first[3])
    evalue = float(first[10])
    return sseqid, pident, alen, evalue


def map_genes_ncbi_rbh(input_csv: str,
                       output_csv: str,
                       email: str,
                       source_taxon: str,
                       target_taxon: str,
                       sleep: float = 0.25,
                       show_progress: bool = True,
                       log_file: str = DEFAULT_LOG) -> Tuple[str, int]:
    """
    NCBI-only RBH mapping flow.

    Parameters:
      - input_csv: CSV with column 'GENE_ID' listing source gene identifiers (one per row)
      - output_csv: path to write mapping results
      - email: contact email used for Entrez/datasets requests
      - source_taxon: free-text taxon name for source proteome (e.g. "Aedes aegypti")
      - target_taxon: free-text taxon name for target proteome (e.g. "Culex pipiens pallens")
      - sleep: pause between per-gene actions
      - show_progress: unused here except for potential future progress hooks
      - log_file: path to append log messages

    Returns (output_csv, number_of_rows)
    """
    logger = _prepare_logger(log_file)
    logger.info("Starting NCBI RBH mapping run")
    logger.info("Input CSV: %s", input_csv)
    logger.info("Source taxon: %s", source_taxon)
    logger.info("Target taxon: %s", target_taxon)

    check_blast_tools()

    # Read input gene IDs
    gene_ids = []
    with open(input_csv) as inf:
        rdr = csv.DictReader(inf)
        for r in rdr:
            gid = r.get("GENE_ID") or r.get("gene_id") or r.get("id")
            gene_ids.append(gid)
    logger.info("Number of input identifiers: %d", len(gene_ids))

    if not NCBI_UTILS_AVAILABLE:
        raise RuntimeError("ncbi_utils.get_proteome_for_taxon not available in repository. Add ncbi_utils.py (datasets/Entrez helper).")

    # Create temporary workdir to store downloads and DBs
    with tempfile.TemporaryDirectory() as workdir:
        logger.info("Using temporary workdir: %s", workdir)

        # Download source proteome
        logger.info("Downloading source proteome for '%s' via NCBI...", source_taxon)
        src_fasta = get_proteome_for_taxon(source_taxon, email, out_dir=workdir)
        if not src_fasta or not os.path.exists(src_fasta):
            raise RuntimeError(f"Failed to download source proteome for '{source_taxon}'")
        logger.info("Source proteome: %s", src_fasta)

        # Download target proteome
        logger.info("Downloading target proteome for '%s' via NCBI...", target_taxon)
        tgt_fasta = get_proteome_for_taxon(target_taxon, email, out_dir=workdir)
        if not tgt_fasta or not os.path.exists(tgt_fasta):
            raise RuntimeError(f"Failed to download target proteome for '{target_taxon}'")
        logger.info("Target proteome: %s", tgt_fasta)

        # Index FASTAs
        logger.info("Indexing source FASTA...")
        src_index = _index_fasta_by_id_and_desc(src_fasta)
        logger.info("Indexing target FASTA...")
        tgt_index = _index_fasta_by_id_and_desc(tgt_fasta)

        # Build BLAST DBs
        logger.info("Building BLAST DB for target FASTA...")
        tgt_db_base = os.path.splitext(tgt_fasta)[0] + "_db"
        make_blast_db(tgt_fasta, tgt_db_base)
        logger.info("Building BLAST DB for source FASTA...")
        src_db_base = os.path.splitext(src_fasta)[0] + "_db"
        make_blast_db(src_fasta, src_db_base)

        # Map provided gene IDs to source SeqRecords
        logger.info("Locating source sequences for provided gene identifiers...")
        gid_to_rec = _find_source_records_for_gene_ids(src_index, gene_ids)

        out_rows = []
        for gid in gene_ids:
            rec = gid_to_rec.get(gid)
            if rec is None:
                logger.warning("No sequence found in source proteome for identifier: %s", gid)
                out_rows.append({
                    "query_id": gid,
                    "source_seq_id": "",
                    "target_seq_id": "",
                    "method_used": "ncbi_rbh",
                    "relationship": "",
                    "percent_identity": "",
                    "coverage": "",
                    "evalue": "",
                    "note": "no_source_sequence_found"
                })
                continue

            # write query seq to temp file
            with tempfile.TemporaryDirectory() as td:
                qpath = os.path.join(td, f"{gid}.fa")
                SeqIO.write(rec, qpath, "fasta")
                blast_out = os.path.join(td, "blast_to_target.tsv")
                logger.info("Blasting source %s -> target DB", rec.id)
                top = _blast_top_hit(qpath, tgt_db_base, blast_out)
                if top is None:
                    logger.info("No hits for %s against target", rec.id)
                    out_rows.append({
                        "query_id": gid,
                        "source_seq_id": rec.id,
                        "target_seq_id": "",
                        "method_used": "ncbi_rbh",
                        "relationship": "no_target_hit",
                        "percent_identity": "",
                        "coverage": "",
                        "evalue": "",
                        "note": "no_hit_to_target"
                    })
                    continue
                tgt_sid, pident, alen, evalue = top
                logger.info("Top target hit for %s -> %s (pident=%.2f evalue=%s)", rec.id, tgt_sid, pident, evalue)

                # find target SeqRecord by id (or substring)
                tgt_rec = tgt_index.get(tgt_sid)
                if tgt_rec is None:
                    # try substring match in headers/descriptions
                    for k, rrec in tgt_index.items():
                        if tgt_sid in k or tgt_sid in (rrec.description or ""):
                            tgt_rec = rrec
                            break

                if tgt_rec is None:
                    logger.warning("Top target hit id %s not found in target FASTA index", tgt_sid)
                    out_rows.append({
                        "query_id": gid,
                        "source_seq_id": rec.id,
                        "target_seq_id": tgt_sid,
                        "method_used": "ncbi_rbh",
                        "relationship": "target_hit_missing_sequence",
                        "percent_identity": pident,
                        "coverage": alen / len(rec.seq) if len(rec.seq) > 0 else "",
                        "evalue": evalue,
                        "note": "target_seq_not_in_fasta_index"
                    })
                    continue

                # now blast the target sequence back against source DB
                t_qpath = os.path.join(td, f"target_{tgt_rec.id}.fa")
                SeqIO.write(tgt_rec, t_qpath, "fasta")
                blast_back_out = os.path.join(td, "blast_back_to_source.tsv")
                logger.info("Blasting target %s -> source DB", tgt_rec.id)
                top_back = _blast_top_hit(t_qpath, src_db_base, blast_back_out)
                if top_back is None:
                    logger.info("No back-hit for %s", tgt_rec.id)
                    out_rows.append({
                        "query_id": gid,
                        "source_seq_id": rec.id,
                        "target_seq_id": tgt_rec.id,
                        "method_used": "ncbi_rbh",
                        "relationship": "no_back_hit",
                        "percent_identity": pident,
                        "coverage": alen / len(rec.seq) if len(rec.seq) > 0 else "",
                        "evalue": evalue,
                        "note": "no_back_hit"
                    })
                    continue
                back_sid, back_pident, back_alen, back_evalue = top_back
                logger.info("Back top hit for %s -> %s", tgt_rec.id, back_sid)

                # Determine if reciprocal
                reciprocal = (back_sid == rec.id) or (rec.id in back_sid) or (back_sid in rec.id)
                if reciprocal:
                    out_rows.append({
                        "query_id": gid,
                        "source_seq_id": rec.id,
                        "target_seq_id": tgt_rec.id,
                        "method_used": "ncbi_rbh",
                        "relationship": "reciprocal_best_hit",
                        "percent_identity": pident,
                        "coverage": alen / len(rec.seq) if len(rec.seq) > 0 else "",
                        "evalue": evalue,
                        "note": ""
                    })
                    logger.info("RBH confirmed: %s <-> %s", rec.id, tgt_rec.id)
                else:
                    out_rows.append({
                        "query_id": gid,
                        "source_seq_id": rec.id,
                        "target_seq_id": tgt_rec.id,
                        "method_used": "ncbi_rbh",
                        "relationship": "best_hit_only",
                        "percent_identity": pident,
                        "coverage": alen / len(rec.seq) if len(rec.seq) > 0 else "",
                        "evalue": evalue,
                        "note": f"back_top={back_sid}"
                    })
            time.sleep(sleep)

        # Write output CSV
        fieldnames = ["query_id", "source_seq_id", "target_seq_id", "method_used", "relationship", "percent_identity", "coverage", "evalue", "note"]
        with open(output_csv, "w", newline="") as outf:
            w = csv.DictWriter(outf, fieldnames=fieldnames)
            w.writeheader()
            for r in out_rows:
                w.writerow(r)
        logger.info("Mapping complete. Wrote %d rows to %s", len(out_rows), output_csv)
        return output_csv, len(out_rows)


# Compatibility wrapper used by the GUI: present so imports succeed.
def map_with_ensembl_and_rbh_ncbi(input_csv: str,
                                  output_csv: str,
                                  email: str,
                                  target_species: str,
                                  auto_download: bool = True,
                                  sleep: float = 0.25,
                                  show_progress: bool = True,
                                  log_file: str = DEFAULT_LOG,
                                  prefer_ncbi: bool = True):
    """
    Backwards-compatible wrapper. In this NCBI-first edition Ensembl-first mapping
    is not implemented — use the NCBI-RBH function directly.

    If this wrapper is called, it raises a clear error explaining how to run NCBI-RBH:
    - call map_genes_ncbi_rbh(input_csv, output_csv, email, source_taxon, target_taxon)
    OR in the GUI enable 'Use NCBI RBH (skip Ensembl)' so the GUI will call map_genes_ncbi_rbh directly.
    """
    raise RuntimeError(
        "map_with_ensembl_and_rbh_ncbi (Ensembl-first wrapper) is not implemented in this NCBI-first build.\n"
        "Please run the NCBI RBH flow instead (map_genes_ncbi_rbh). Example:\n"
        "  from ensembl_rbh import map_genes_ncbi_rbh\n"
        "  map_genes_ncbi_rbh('input.csv','output.csv','you@example.com','Aedes aegypti','Culex pipiens pallens')\n"
        "Or enable 'Use NCBI RBH' in the GUI so it will call the NCBI RBH path."
    )

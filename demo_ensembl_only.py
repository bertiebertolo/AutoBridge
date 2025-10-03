#!/usr/bin/env python3
"""
Quick demo: query Ensembl homology for a small set of AAEL gene IDs
Writes results to demo_ensembl_results.csv

No BLAST required. This demonstrates the Ensembl lookup part of AUTOBridge.
"""
import csv
import requests
import time

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Accept": "application/json", "User-Agent": "AUTOBridge-demo your.email@example.com"}

# small demo input (replace or load from a CSV)
demo_input = [
    {"GENE_ID": "AAEL000020", "HD": "FALSE", "JO": "TRUE"},
    {"GENE_ID": "AAEL000076", "HD": "FALSE", "JO": "TRUE"},
    {"GENE_ID": "AAEL000084", "HD": "FALSE", "JO": "TRUE"},
]

TARGET_SPECIES_SUBSTR = "culex"  # change to any species string you want

def query_homology(gene_id, target_species_substr):
    url = f"{ENSEMBL_REST}/homology/id/{gene_id}"
    params = {"type": "orthologues", "sequence": 0, "format": "full"}
    r = requests.get(url, params=params, headers=HEADERS, timeout=30)
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
                "query_id": gene_id,
                "target_species": species,
                "target_gene_id": tgt.get("gene_id") or tgt.get("id",""),
                "target_protein_id": tgt.get("id",""),
                "orthology_type": h.get("type",""),
                "confidence": h.get("confidence",""),
                "perc_id": h.get("percentage_identities") or h.get("perc_id","")
            })
    return out

def main():
    out_rows = []
    for rec in demo_input:
        gid = rec["GENE_ID"]
        try:
            homs = query_homology(gid, TARGET_SPECIES_SUBSTR)
        except Exception as e:
            out_rows.append({
                "query_id": gid,
                "target_gene_id": "",
                "target_species": "",
                "orthology_type": "",
                "confidence": "",
                "perc_id": "",
                "HD": rec["HD"],
                "JO": rec["JO"],
                "note": f"ensembl_error: {e}"
            })
            continue

        if not homs:
            out_rows.append({
                "query_id": gid,
                "target_gene_id": "",
                "target_species": "",
                "orthology_type": "",
                "confidence": "",
                "perc_id": "",
                "HD": rec["HD"],
                "JO": rec["JO"],
                "note": "no_ensembl_hit"
            })
        else:
            # write one row per homology returned (could be many)
            for h in homs:
                out_rows.append({
                    "query_id": h["query_id"],
                    "target_gene_id": h["target_gene_id"],
                    "target_species": h["target_species"],
                    "orthology_type": h["orthology_type"],
                    "confidence": h["confidence"],
                    "perc_id": h["perc_id"],
                    "HD": rec["HD"],
                    "JO": rec["JO"],
                    "note": "ensembl"
                })
        # be polite
        time.sleep(0.2)

    # write CSV
    out_csv = "demo_ensembl_results.csv"
    fieldnames = ["query_id","target_gene_id","target_species","orthology_type","confidence","perc_id","HD","JO","note"]
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in out_rows:
            writer.writerow(r)
    print(f"Wrote {len(out_rows)} rows to {out_csv}")

if __name__ == "__main__":
    main()
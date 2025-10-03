#!/usr/bin/env python3
"""
Map Aedes aegypti AAEL gene IDs to Culex orthologs (Ensembl REST) with file pickers.

Usage:
    python ensembl_map_aael_to_culex_gui.py

What it does:
- Opens a file dialog to select the input CSV (expects columns: GENE_ID, HD, JO).
- Asks for contact email and target-species substring (e.g. "culex" or "culex_pipiens_pallens").
- Opens a Save dialog to choose output CSV path/name.
- Queries Ensembl REST /homology/id/<gene> for orthologues and writes mappings to the chosen CSV.

Dependencies:
    pip install pandas requests
Tkinter is part of the Python standard library (usually installed by default).

Notes:
- The script prefers one-to-one orthologs by default; you can change that in the dialog.
- For maximum accuracy you can set a longer sleep between requests.
"""

import sys
import time
import csv
import requests
import pandas as pd
from typing import List, Dict
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox

ENSEMBL_REST = "https://rest.ensembl.org"

HEADERS_TEMPLATE = {"Accept": "application/json"}

def make_session(retries=5, backoff_factor=0.5, status_forcelist=(500,502,503,504)):
    session = requests.Session()
    retry = Retry(total=retries, read=retries, connect=retries,
                  backoff_factor=backoff_factor,
                  status_forcelist=status_forcelist,
                  allowed_methods=frozenset(['GET','POST']))
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

def query_homology(session: requests.Session, gene_id: str, user_agent: str,
                   target_species_substr: str = None, timeout=30) -> List[Dict]:
    headers = HEADERS_TEMPLATE.copy()
    headers["User-Agent"] = user_agent
    url = f"{ENSEMBL_REST}/homology/id/{gene_id}"
    params = {"type": "orthologues", "sequence": 0, "format": "full"}
    resp = session.get(url, headers=headers, params=params, timeout=timeout)
    if resp.status_code == 429:
        raise requests.HTTPError("429 rate limit")
    resp.raise_for_status()
    j = resp.json()
    homologies_out = []
    for rec in j.get("data", []):
        for h in rec.get("homologies", []):
            tgt = h.get("target", {})
            species = tgt.get("species", "")
            if target_species_substr:
                if target_species_substr.lower() not in species.lower():
                    continue
            # prefer gene_id when available
            target_gene_id = tgt.get("gene_id") or tgt.get("id") or ""
            target_protein_id = tgt.get("id") if tgt.get("id") != target_gene_id else ""
            homologies_out.append({
                "target_gene_id": target_gene_id,
                "target_protein_id": target_protein_id,
                "target_species": species,
                "orthology_type": h.get("type", ""),
                "confidence": h.get("confidence", ""),
                "perc_id": h.get("percentage_identities") or h.get("percent_id") or ""
            })
    return homologies_out

def map_genes(input_path: str, output_path: str, target_species: str, email: str,
              sleep: float = 0.25, prefer_one2one: bool = True):
    df = pd.read_csv(input_path)
    if "GENE_ID" not in df.columns:
        raise ValueError("Input CSV must have a 'GENE_ID' column")
    session = make_session()
    user_agent = f"ensembl-orthology-script ({email})" if email else "ensembl-orthology-script"
    results = []
    total = len(df)
    for idx, row in df.iterrows():
        query_id = str(row["GENE_ID"]).strip().strip('"')
        hd = row.get("HD", "")
        jo = row.get("JO", "")
        backoff = 1.0
        max_backoff = 60.0
        got = None
        while True:
            try:
                homs = query_homology(session, query_id, user_agent, target_species_substr=target_species)
                got = homs
                break
            except requests.HTTPError as he:
                print(f"[{idx+1}/{total}] HTTP error for {query_id}: {he}. Backing off {backoff}s...", file=sys.stderr)
                time.sleep(backoff)
                backoff = min(max_backoff, backoff * 2)
            except Exception as e:
                print(f"[{idx+1}/{total}] Error querying {query_id}: {e}", file=sys.stderr)
                got = [{"error": str(e)}]
                break

        if not got:
            results.append({
                "query_id": query_id,
                "target_gene_id": "",
                "target_protein_id": "",
                "target_species": "",
                "orthology_type": "",
                "confidence": "",
                "perc_id": "",
                "HD": hd,
                "JO": jo,
                "note": "no_homology_found"
            })
        else:
            if len(got) == 1 and got[0].get("error"):
                results.append({
                    "query_id": query_id,
                    "target_gene_id": "",
                    "target_protein_id": "",
                    "target_species": "",
                    "orthology_type": "",
                    "confidence": "",
                    "perc_id": "",
                    "HD": hd,
                    "JO": jo,
                    "note": got[0].get("error")
                })
            else:
                selected = []
                if prefer_one2one:
                    for h in got:
                        typ = h.get("orthology_type","").lower()
                        if "one2one" in typ or "one_to_one" in typ or "one-to-one" in typ:
                            selected.append(h)
                if not selected:
                    selected = got
                for h in selected:
                    results.append({
                        "query_id": query_id,
                        "target_gene_id": h.get("target_gene_id", ""),
                        "target_protein_id": h.get("target_protein_id", ""),
                        "target_species": h.get("target_species", ""),
                        "orthology_type": h.get("orthology_type", ""),
                        "confidence": h.get("confidence", ""),
                        "perc_id": h.get("perc_id", ""),
                        "HD": hd,
                        "JO": jo,
                        "note": ""
                    })
        time.sleep(sleep)

    out_fields = ["query_id","target_gene_id","target_protein_id","target_species",
                  "orthology_type","confidence","perc_id","HD","JO","note"]
    with open(output_path, "w", newline="") as outf:
        writer = csv.DictWriter(outf, fieldnames=out_fields)
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    return len(results)

def run_gui():
    root = tk.Tk()
    root.withdraw()
    messagebox.showinfo("Info", "You will be asked to select the input CSV, then choose where to save the output CSV.")
    input_path = filedialog.askopenfilename(title="Select input CSV", filetypes=[("CSV files","*.csv"),("All files","*.*")])
    if not input_path:
        messagebox.showwarning("Canceled", "No input file selected. Exiting.")
        return
    save_path = filedialog.asksaveasfilename(title="Save output CSV as", defaultextension=".csv",
                                             filetypes=[("CSV files","*.csv"),("All files","*.*")],
                                             initialfile="aael_to_culex_mapping.csv")
    if not save_path:
        messagebox.showwarning("Canceled", "No output file selected. Exiting.")
        return
    # required email
    email = None
    while not email:
        email = simpledialog.askstring("Contact email", "Enter contact email (required for User-Agent):")
        if email is None:
            messagebox.showwarning("Canceled", "Email is required. Exiting.")
            return
        email = email.strip()
    target_species = simpledialog.askstring("Target species", "Enter substring to match target species (default 'culex'):", initialvalue="culex")
    if target_species is None:
        target_species = "culex"
    try:
        sleep = simpledialog.askfloat("Sleep (seconds)", "Seconds to sleep between requests (default 0.25):", initialvalue=0.25)
        if sleep is None:
            sleep = 0.25
    except Exception:
        sleep = 0.25
    pref = messagebox.askyesno("Prefer one-to-one?", "Prefer one-to-one orthologs when available? (Yes = prefer one-to-one)", default=messagebox.YES)
    # Run mapping
    try:
        nrows = map_genes(input_path, save_path, target_species, email, sleep=sleep, prefer_one2one=pref)
        messagebox.showinfo("Done", f"Mapping complete. Wrote {nrows} row(s) to:\n{save_path}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred:\n{e}")

if __name__ == "__main__":
    run_gui()
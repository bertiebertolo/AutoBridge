### AUTOBridge

AUTOBridge is a lightweight desktop GUI and CLI toolkit to map orthologs and run BLAST/RBH workflows across species.  
It supports two primary data sources for proteomes and orthology information:

- Ensembl (preferred for curated homology calls where available)
- NCBI (datasets / Taxonomy / Assembly FTP) — used for proteome downloads and NCBI-only RBH mapping

Recommendation: use both resources where possible. Ensembl provides curated orthology for many model organisms; NCBI provides broader assembly coverage (including strains/subspecies) and is used for automatic proteome download + RBH fallback when Ensembl lacks the species or homologies. 


---

## Highlights (what AUTOBridge does)
- Prefer Ensembl homology lookups (fast, curated one‑to‑one orthologues).
- RBH fallback using BLAST+ and a target proteome FASTA when Ensembl lacks a species or a homology.
- NCBI-first RBH mode: automatically download proteomes from NCBI using the Datasets CLI (recommended) or Entrez fallback, then run Reciprocal Best Hit (RBH) mapping.
- GUI (Tkinter) that:
  - suggests species names (NCBI Taxonomy first, Ensembl fallback),
  - tails autobridge.log for live progress,
  - exposes options: Ensembl-first, NCBI auto-download, and NCBI-only RBH.
- Outputs reproducible CSV mappings with provenance and BLAST/orthology metrics.

---

## Table of contents
1. Requirements
2. Installation (recommended: conda)
3. Manual / alternative installs
4. Verify your environment
5. Quick start — GUI
6. Quick start — command line
7. Input / output formats
8. How orthology is determined (Ensembl + RBH)
9. Choosing species & why use both NCBI and Ensembl
10. Tips & best practices
11. Troubleshooting
12. Files in the repository
13. Contributing & license

---

## 1) Requirements
- Python 3.8+ (3.10 recommended)
- BLAST+ executables on PATH: makeblastdb, blastp, tblastn, blastn
- Network access for Ensembl/NCBI calls and downloads
- Recommended (not strictly required):
  - NCBI Datasets CLI (ncbi-datasets-cli) — preferred for robust proteome downloads
  - Biopython — Entrez fallback and FASTA parsing
  - requests, tqdm (Python libraries)

---

## 2) Installation (recommended: conda)
A reproducible conda environment is recommended (bioconda & conda-forge):

```bash
# Create environment and add channels (if required)
conda create -n autobridge -y python=3.10
conda activate autobridge
conda config --add channels conda-forge
conda config --add channels bioconda

# Install executables and Python packages
conda install -c bioconda blast ncbi-datasets-cli -y
conda install -c conda-forge biopython requests pandas tqdm -y
```

Notes:
- `blast` (BLAST+) provides makeblastdb, blastp, tblastn, blastn.
- `ncbi-datasets-cli` is strongly recommended; if absent, AUTOBridge will try an Entrez/FTP fallback but results may be less reliable.

---

## 3) Manual / alternative installs
- macOS Homebrew:
  ```bash
  brew install ncbi-datasets-cli
  brew install blast
  pip install biopython requests pandas tqdm
  ```
- Official NCBI BLAST binaries:
  - Download from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  - Extract and add `bin/` to PATH.
- If you cannot install `ncbi-datasets-cli`, ensure Biopython is installed and provide a valid Entrez email in the GUI.

---

## 4) Verify your environment
Run these commands to confirm essentials:

```bash
which datasets && datasets --version
which makeblastdb && makeblastdb -version
python -c "import Bio,requests,pandas; print('python deps OK')"
python -c "import tkinter; print('tkinter OK')"
```

If any command fails, address that dependency before running mapping.

---

## 5) Quick start — GUI
1. Launch:
   ```bash
   python autobridge_gui.py
   ```
2. In the Ortholog Mapping tab:
   - Input CSV: select CSV with a `GENE_ID` column (one identifier per row).
   - Output CSV: choose an output path for results.
   - Source taxon: the species for your input IDs (e.g., `Aedes aegypti`).
   - Target taxon: the species you want orthologues in (e.g., `Culex pipiens pallens`).
   - Entrez email: set your contact email (used in Entrez and User-Agent strings).
   - Options:
     - "Auto-download from NCBI": allow automatic NCBI proteome downloads if needed.
     - "Use NCBI RBH (skip Ensembl)": force NCBI-only RBH mapping (useful when Ensembl lacks the target).
     - "Show progress output": show tqdm/simple progress.
   - Click "Run Mapping". The GUI will tail `autobridge.log` so you can watch download, makeblastdb and BLAST progress.

---

## 6) Quick start — command line

- Ensembl-first (if Ensembl supports both species):
```bash
python - <<'PY'
from ensembl_rbh import map_genes
out,n = map_genes("input.csv","output.csv","you@example.com","culex_quinquefasciatus", show_progress=True)
print(out,n)
PY
```

- NCBI-only RBH (downloads proteomes and runs RBH):
```bash
python - <<'PY'
from ensembl_rbh import map_genes_ncbi_rbh
out,n = map_genes_ncbi_rbh("input.csv","output.csv","you@example.com","Aedes aegypti","Culex pipiens pallens")
print(out,n)
PY
```

Notes:
- Always supply a valid Entrez email for polite requests to NCBI and Ensembl.
- For reproducibility, you can manually download FASTA files and adapt the mapping calls to pass local FASTA paths.

---

## 7) Input / output formats

Input CSV:
- Required column: `GENE_ID` (case-insensitive `gene_id` or `id` accepted).
- Optional columns `HD`, `JO` are preserved in outputs where relevant.

Example:
```
GENE_ID,HD,JO
AAEL000020,TRUE,FALSE
AAEL000100,TRUE,FALSE
```

Output CSV (NCBI-RBH / combined):
- query_id
- source_seq_id
- target_seq_id
- method_used (ensembl | ncbi_rbh | rbh_low_quality | ensembl_no_rbh, ...)
- relationship (reciprocal_best_hit, best_hit_only, no_target_hit, no_back_hit, ...)
- percent_identity
- coverage
- evalue
- note

---

## 8) How orthology is determined (Ensembl + RBH)
Default behaviour (per gene):

1. Ensembl-first:
   - Query Ensembl REST /homology for the input gene.
   - If curated orthologue(s) exist (prefer one-to-one), record them with method_used=`ensembl`.

2. RBH fallback (if Ensembl lacks homology or if NCBI-only mode selected):
   - Obtain the query protein sequence:
     - Prefer local source proteome FASTA (if provided),
     - Otherwise fetch sequence from Ensembl sequence endpoint or from downloaded source proteome.
   - BLAST (blastp) the query against the target proteome DB.
   - If the top hit passes thresholds, retrieve the target sequence and BLAST it back against the source DB.
   - If reciprocal top hit returns the original query sequence (or matches closely), call it `reciprocal_best_hit`.

Default thresholds (configurable in code):
- evalue <= 1e-5
- percent identity >= 30%
- coverage >= 0.5 (50% of query)

Caveats:
- RBH is heuristic, not a formal phylogenetic orthology call; curate results for paralogous families.
- If your input IDs are gene symbols (not present in FASTA headers), a symbol→protein ID mapping step may be required.

---

## 9) Choosing species & why use both NCBI and Ensembl
- Ensembl: curated orthology for many model organisms; use canonical Ensembl names (e.g. `aedes_aegypti`, `culex_quinquefasciatus`) when relying on Ensembl REST.
- NCBI: broader assembly coverage (including subspecies/strains). Use free-text taxon names (e.g., `Culex pipiens pallens`) and enable Auto-download to fetch proteomes via NCBI Datasets.

Recommendation:
- Try Ensembl-first for genes/species covered by Ensembl. When Ensembl lacks coverage or specific assemblies, enable NCBI Auto-download or run NCBI RBH. Using both improves coverage and reduces missed orthologues.

Suggested workflow:
1. Run Ensembl homology for your list (fast and curated).
2. For genes missing Ensembl homology, run NCBI RBH against a downloaded proteome for your target species.

---

## 10) Tips & best practices
- Pre-download/store proteome FASTAs and pre-build BLAST DBs if you run repeated analyses.
- For large jobs, run Ensembl homology in batch first; only RBH the unresolved subset.
- Use multiple threads in BLAST for speed (modify blast_utils if desired).
- Avoid spaces in file paths.
- To debug downloads, set `ncbi_utils` to write to a persistent folder and inspect the extracted contents.

---

## 11) Troubleshooting
- "Failed to download proteome": ensure `datasets` is installed and on PATH (`which datasets`). If unavailable, the Entrez fallback may still work if Biopython is installed and a valid Entrez email is provided.
- BLAST not found: install BLAST+ and ensure `makeblastdb` and `blastp` are on PATH.
- Ensembl / NCBI rate-limits: include an email and respect polite sleeps (the code includes sleep between requests).
- If network/proxy blocks downloads: set HTTP_PROXY/HTTPS_PROXY environment variables or run from a network with access.

---

## 12) Files in the repository (where to edit)
- `autobridge_gui.py` — Tkinter GUI main
- `ensembl_rbh.py` — mapping logic; includes `map_genes` (Ensembl-first) and `map_genes_ncbi_rbh` (NCBI-only RBH)
- `ncbi_utils.py` — NCBI download helpers (datasets CLI preferred; Entrez fallback)
- `blast_utils.py` — BLAST wrappers (make DB, call blastp/tblastn)
- `debug_demo.py` — small helper for connectivity/BLAST checks
- `requirements.txt` — Python dependency list (if provided)

Edit these files to adjust RBH thresholds, logging, or to add persistent download options.

---

## 13) Contributing & license
- License: MIT-style by default (add `LICENSE` to the repo to make it explicit).
- Contributions: open issues or PRs. If you want, maintainers can prepare a PR branch `feature/ncbi-first` containing NCBI-first changes.

---

```

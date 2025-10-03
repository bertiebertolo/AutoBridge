
### AUTOBridge

AUTOBridge is a desktop GUI and CLI toolkit to automate ortholog discovery and BLAST workflows across species. It queries Ensembl/VectorBase for homologies (preferred), and when Ensembl does not provide a clear one‑to‑one ortholog it falls back to a reciprocal-best-hit (RBH) BLAST procedure for higher confidence. AUTOBridge also provides ad-hoc protein and gene BLAST tabs and supports optional automatic proteome download from NCBI (RefSeq preferred).

Short description
AUTOBridge — GUI to automatically map orthologs and run BLAST (Ensembl + RBH fallback)

Highlights
- Ensembl homology lookup (fast; prefers curated one‑to‑one orthologs)
- RBH fallback using BLAST+ for missing/ambiguous cases
- Protein BLAST (blastp) and Gene BLAST (blastn / tblastn / blastx) GUI tabs
- Optional automatic proteome download from NCBI (Entrez)
- Outputs reproducible CSV mappings with provenance and BLAST metrics
- Minimal dependencies; conda recommended for reproducible installs

Quick index
- Installation (conda) — recommended
- Manual NCBI BLAST+ placement — optional
- Basic usage (GUI & debug helper)
- Gene BLAST vs Protein BLAST — which command to use
- RBH / ortholog mapping details
- Tips & troubleshooting

---

## 1) Install (recommended: conda)

Recommended single-line (creates env and installs Python packages + BLAST+):
```bash
conda create -n autobridge -c conda-forge -c bioconda python=3.10 pandas requests biopython blast -y
conda activate autobridge
```

Step-by-step:
1. Create & activate environment:
   ```bash
   conda create -n autobridge python=3.10 -y
   conda activate autobridge
   ```
2. Install Python packages:
   ```bash
   conda install -c conda-forge pandas requests biopython -y
   ```
3. Install BLAST+ (bioconda):
   ```bash
   conda install -c bioconda blast -y
   ```

If you prefer pip inside conda:
```bash
conda create -n autobridge python=3.10 -y
conda activate autobridge
conda install -c conda-forge pip -y
pip install -r requirements.txt
conda install -c bioconda blast -y
```

Verify:
```bash
python -c "import pandas,requests,Bio; print('python deps OK')"
which makeblastdb && which blastp && which blastdbcmd
python -c "import tkinter; print('tkinter OK')"   # GUI check
```

---

## 2) Manual BLAST+ installation (optional)

If you prefer the official NCBI binary distribution:
1. Download the BLAST+ toolkit (2.16.0+ recommended) from:
   https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
2. Extract — folder will be named like `ncbi-blast-2.16.0+`.
3. Place the extracted `ncbi-blast-2.16.0+` folder next to your project (e.g., inside `AUTOBridge/`) or anywhere you like.

Two ways to use it:
- Call full path to executables:
  ```bash
  ./AUTOBridge/ncbi-blast-2.16.0+/bin/blastp -version
  ```
- Or add the `bin/` folder to your PATH (recommended):
  - macOS / Linux:
    ```bash
    export PATH="/full/path/to/ncbi-blast-2.16.0+/bin:$PATH"
    ```
    Add that line to `~/.bashrc`/`~/.zshrc` to persist.
  - Windows (PowerShell temporary):
    ```powershell
    $env:Path += ";C:\full\path\to\ncbi-blast-2.16.0+\bin"
    ```
    Or add via System > Environment Variables.

If a binary is flagged as untrusted you may need to unblock it:
- macOS: `xattr -d com.apple.quarantine /path/to/bin/blastp`
- Linux: `chmod +x /path/to/bin/*`
- Windows: Properties → Unblock or run as Admin

Note: using conda/bioconda is simpler because BLAST becomes available automatically in the conda env PATH.

---

## 3) Basic usage

Start the GUI:
```bash
python autobridge_gui.py
```

Quick debug helper (small check of connectivity and BLAST presence):
```bash
python debug_demo.py --gene AAEL000020 --species culex --email you@example.com
```

Input format for ortholog mapping:
- CSV with header including at least: `GENE_ID` (AAEL IDs), optional `HD`, `JO`
- Example row: `AAEL000020,TRUE,FALSE`

Output:
- CSV with columns:
  - query_id (AAEL)
  - target_gene_id (Culex / target species ID found)
  - method_used (ensembl | rbh | rbh_failed | no_rbh_available | etc.)
  - orthology_type (from Ensembl if available)
  - percent_identity (from BLAST if used)
  - coverage (alignment coverage fraction)
  - evalue (BLAST evalue)
  - HD, JO (copied)
  - note (additional info/errors)

---

## 4) Gene BLAST vs Protein BLAST — which to use?

Yes — BLAST supports both gene (nucleotide) and protein searches. Use the program that matches your data and goals:

- blastp (protein → protein)
  - Use when both query and subject are protein sequences (FASTA of amino acids).
  - Typical for RBH orthology at protein level.

- blastn (nucleotide → nucleotide)
  - Use when both query and subject are nucleotide sequences (DNA/RNA).
  - Works for comparing gene sequences directly.

- tblastn (protein query → translated nucleotide DB)
  - Use when you have protein queries and the subject DB is nucleotide (genome/transcriptome). BLAST translates the DB in six frames and searches proteins.
  - Useful to find genes in a genome using a protein query.

- blastx (translated query → protein DB)
  - Use when query is nucleotide (e.g., transcript) and you want to search a protein DB (the query is translated in six frames).

- tblastx (translated query → translated DB)
  - Translates both query and DB and compares proteins; heavier and slower, rarely required.

AUTOBridge GUI exposes:
- Protein BLAST tab (blastp)
- Gene BLAST tab (blastn / tblastn style — depending on whether subject DB is nucleotide or protein you may choose tblastn)

Example commands:
```bash
# make protein DB
makeblastdb -in culex_proteins.faa -dbtype prot -out culex_db

# run blastp (protein vs protein)
blastp -query query.faa -db culex_db -outfmt 6 -max_target_seqs 5 -evalue 1e-5 -out blastp_out.tsv

# make nucleotide DB
makeblastdb -in genome.fa -dbtype nucl -out genome_db

# run tblastn (protein query vs nucleotide DB translated)
tblastn -query query.faa -db genome_db -outfmt 6 -max_target_seqs 5 -evalue 1e-5 -out tblastn_out.tsv

# blastn (nuc vs nuc)
blastn -query query.fa -db subject_nuc_db -outfmt 6 -out blastn_out.tsv
```

---

## 5) RBH ortholog mapping (how AUTOBridge determines orthologs)

AUTOBridge mapping flow (per input gene):
1. Try Ensembl homology endpoint (/homology/id/<gene>) — accept curated one-to-one orthologs when available.
2. If Ensembl returns no clear ortholog:
   - Obtain the protein sequence for the AAEL gene:
     - Prefer a local Aedes proteome FASTA (fast).
     - Or fetch from Ensembl REST (slower).
   - BLAST the Aedes protein against the target species proteome (blastp).
   - If top hit passes quality thresholds (defaults: evalue ≤ 1e-5, pct identity ≥ 30%, coverage ≥ 50%), retrieve the subject sequence.
   - BLAST the subject back against the Aedes proteome (reciprocal blast).
   - If the reciprocal top hit maps back to the original AAEL ID (or matches closely), call it an RBH and record as an ortholog.
   - Otherwise mark as ambiguous (record candidates and metrics).

Default RBH thresholds (configurable in GUI):
- evalue: 1e-5
- pct identity: 30%
- coverage: 50% (fraction of query)

Notes: RBH is a strong heuristic but not a formal phylogenetic orthology inference — combining Ensembl curated calls with RBH gives practical coverage and accuracy.

---

## 6) Choosing target species in the GUI

- Preferred: use canonical Ensembl species names (lowercase, underscores), e.g.:
  - `homo_sapiens`, `aedes_aegypti`, `anopheles_gambiae`, `culex_pipiens`
- Flexible: you can enter a simple substring (e.g., `culex`) and AUTOBridge will substring-match Ensembl species entries.
- Use the "Suggest species" button in the GUI to query Ensembl's `/info/species` and get canonical names to paste into the species field.

---

## 7) Tips & best practices

- Prebuild BLAST DBs for target proteomes and reuse them:
  ```bash
  makeblastdb -in culex_proteins.faa -dbtype prot -out culex_db
  ```
- For thousands of genes:
  - Run Ensembl homology in batch first (fast).
  - Only run RBH for the subset missing clear Ensembl orthologs.
- Use `-num_threads N` for faster BLAST on multi-core systems (you can add this in code or modify GUI to expose it).
- Avoid spaces in BLAST filenames/paths.
- If many paralogs exist, RBH may return paralogous matches — consider further filtering or manual curation for gene families.

---

## 8) Troubleshooting

- If `blastp` or `makeblastdb` not found:
  - Ensure conda env is activated or the NCBI BLAST bin/ is in PATH.
  - On macOS/Linux check `which blastp`. On Windows use PowerShell `Get-Command blastp`.
- BLAST flagged as untrusted:
  - Unblock on Windows; use `chmod +x` on Linux/macOS; on macOS remove quarantine using `xattr -d`.
- Ensembl rate limits / connectivity:
  - Include a contact email in User-Agent when calling Ensembl REST.
  - The scripts include polite sleeps and retry/backoff; you can increase the sleep in GUI settings.
- If you run into package resolution issues in conda, use `mamba` as a faster solver:
  ```bash
  conda install -n base -c conda-forge mamba
  mamba create -n autobridge -c conda-forge -c bioconda python=3.10 pandas requests biopython blast -y
  ```

---

## 9) Included helper scripts & files
- `autobridge_gui.py` — main GUI (Tkinter)
- `ensembl_rbh.py` — Ensembl + RBH logic
- `blast_utils.py` — BLAST wrappers
- `ncbi_utils.py` — optional automatic NCBI proteome download helpers
- `debug_demo.py` — minimal environment & connectivity check
- `requirements.txt` — Python package list (pandas, requests, biopython)

---

## 10) License & contribution
Add your preferred license and contribution guidelines in `LICENSE` and `CONTRIBUTING.md` as desired. If you want, I can prepare default MIT or Apache-2.0 license files.

---

If you want, I can:
- produce a one-line README short description optimized for GitHub's repo description field,
- add an example input CSV and a small test dataset for quick verification,
- or prepare a headless CLI wrapper for batch runs on HPC (non-GUI).

```

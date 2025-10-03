```markdown
# AUTOBridge

AUTOBridge is a desktop GUI application to automatically map genes across species (ortholog mapping) and run gene/protein BLASTs.

New: specify any target species by typing it in the GUI and use "Suggest species" to query Ensembl's species index.

Quick guidelines for the "target species" field
- Recommended format: canonical Ensembl species name (lowercase, underscore-separated), e.g.:
  - Homo sapiens -> `homo_sapiens`
  - Aedes aegypti -> `aedes_aegypti`
  - Anopheles gambiae -> `anopheles_gambiae`
  - Culex pipiens -> `culex_pipiens`
  - Culex pipiens pallens -> `culex_pipiens_pallens` (if present)
- Flexible option: enter a substring (case-insensitive), e.g. `culex` or `aegypti`. The tool will match substrings.
- Use "Suggest species" in the GUI to query Ensembl's species list and copy the exact Ensembl `name` if you prefer the canonical form.

Conda-based install (recommended)
- Single one-line (preferred if you want everything in one step):
  - conda create -n autobridge -c conda-forge -c bioconda python=3.10 pandas requests biopython blast -y
  - conda activate autobridge
- Step-by-step (clearer if you want more control):
  1. Create and activate environment:
     - conda create -n autobridge python=3.10 -y
     - conda activate autobridge
  2. Install Python packages:
     - conda install -c conda-forge pandas requests biopython -y
  3. Install BLAST+ (from bioconda):
     - conda install -c bioconda blast -y

If you prefer to use pip inside the conda env:
1. conda create -n autobridge python=3.10 -y
2. conda activate autobridge
3. conda install -c conda-forge pip -y
4. pip install -r requirements.txt
5. conda install -c bioconda blast -y

Verify install
- Python deps:
  python -c "import pandas,requests,Bio; print('python deps OK')"
- BLAST tools:
  which makeblastdb && which blastp && which blastdbcmd
  (or on Windows PowerShell: Get-Command makeblastdb, Get-Command blastp)

BLAST notes & best practices
- For RBH (reciprocal BLAST) accuracy, provide local proteome FASTAs (Culex proteome and Aedes proteome).
- Build BLAST DBs once with makeblastdb and reuse them:
  makeblastdb -in culex_proteins.faa -dbtype prot -out culex_db
- For speed, add `-num_threads N` to blastp when running on multi-core machines.

Debug helper (quick)
- Use the provided `debug_demo.py` to check:
  - BLAST+ availability,
  - Ensembl REST /info/species connectivity,
  - A single homology lookup for one AAEL ID.
- Example:
  python debug_demo.py --gene AAEL000020 --species culex --email you@example.com

Removed demo
- I removed the longer demo script from the repository and replaced it with the minimal `debug_demo.py` above. The rest of AUTOBridge (GUI and core modules) is unchanged.

Running AUTOBridge GUI
- Install conda env as above and ensure BLAST is on PATH.
- Run:
  python autobridge_gui.py

If you'd like a headless CLI runner or packaged executables (Windows .exe / macOS app), I can produce those next.
```
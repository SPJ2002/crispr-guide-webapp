# CRISPR Guide Design Web App

A small full-stack tool that finds SpCas9 guide RNAs in a DNA sequence and ranks them by simple quality and off-target checks.

## Features
- Paste a DNA sequence (A/C/G/T, optional FASTA header, up to 10,000 bp).
- Finds every 20-nt guide followed by an NGG PAM on **both strands**.
- Reports GC content, flags (poly-T runs, GC outside 40-70%) and a low/medium/high risk label.
- Counts **off-target near-matches**: other 20-mers in the sequence (either strand) within 3 mismatches. Implemented with one-hot encoding and matrix multiplication (NumPy), so a 10 kb sequence with about 1,300 guides runs in under 0.2 s.
- Stores each run in SQLite and shows recent runs; the results table is sortable.

## Stack
Python, Flask (REST API), SQLite, NumPy, vanilla JavaScript/HTML/CSS, pytest.

## API
- `POST /api/analyze` with JSON `{"sequence": "..."}` returns guides and summary (400 for invalid input, 413 if too long)
- `GET /api/history` returns the last 10 runs

## Run
```bash
pip install -r requirements.txt
python app.py        # open http://127.0.0.1:5000
python -m pytest -q  # 13 tests
```

## Limitations
- Off-target search only looks inside the pasted sequence, not a whole genome, and ignores PAM requirements when counting near-matches. Real design tools (e.g. CRISPOR, Cas-OFFinder) search a full reference genome.
- The risk label is a simple heuristic, not a validated efficiency score.
- Next steps: search against a reference genome index, add on-target scoring, and deploy the app.

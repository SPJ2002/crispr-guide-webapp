"""Core logic: guide finding, GC/flag checks and mismatch-based off-target search."""
import re
import numpy as np

COMP = str.maketrans("ACGT", "TGCA")
GUIDE_LEN = 20
MAX_MISMATCH = 3


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def clean(raw):
    """Strip FASTA headers/whitespace, uppercase, and validate A/C/G/T only."""
    lines = [l for l in raw.splitlines() if not l.startswith(">")]
    seq = re.sub(r"\s", "", "".join(lines)).upper()
    if not seq or not re.fullmatch(r"[ACGT]+", seq):
        raise ValueError("Sequence must contain only A, C, G, T")
    return seq


def gc_content(guide):
    return round((guide.count("G") + guide.count("C")) / len(guide), 3)


def find_guides(seq):
    """Find 20-nt guides followed by an NGG PAM (SpCas9) on both strands."""
    guides = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for i in range(len(s) - GUIDE_LEN - 2):
            if s[i + GUIDE_LEN + 1:i + GUIDE_LEN + 3] == "GG":
                start = i if strand == "+" else len(seq) - (i + GUIDE_LEN)
                guides.append({"guide": s[i:i + GUIDE_LEN], "pam": s[i + GUIDE_LEN:i + GUIDE_LEN + 3],
                               "strand": strand, "start": start})
    return guides


def _onehot(arr):
    out = np.zeros((len(arr), GUIDE_LEN * 4), dtype=np.float32)
    for k in range(GUIDE_LEN):
        out[np.arange(len(arr)), k * 4 + arr[:, k]] = 1
    return out


def count_offtargets(seq, guides, batch=500):
    """Count other 20-mers (both strands) within MAX_MISMATCH mismatches of each guide.
    Uses one-hot encoding and matrix multiplication so it stays fast."""
    if len(seq) < GUIDE_LEN or not guides:
        return [0] * len(guides)
    code = np.frombuffer(seq.encode(), dtype=np.uint8)
    lut = np.zeros(256, dtype=np.int64)
    for i, b in enumerate(b"ACGT"):
        lut[b] = i
    fwd = np.lib.stride_tricks.sliding_window_view(lut[code], GUIDE_LEN)
    rc = np.lib.stride_tricks.sliding_window_view(lut[np.frombuffer(revcomp(seq).encode(), dtype=np.uint8)], GUIDE_LEN)
    windows = _onehot(np.vstack([fwd, rc]))
    garr = np.array([[lut[ord(c)] for c in g["guide"]] for g in guides])
    counts = []
    for s in range(0, len(garr), batch):
        matches = windows @ _onehot(garr[s:s + batch]).T
        counts.extend(((GUIDE_LEN - matches) <= MAX_MISMATCH).sum(axis=0) - 1)  # minus the guide itself
    return [int(max(c, 0)) for c in counts]


def analyze(raw):
    seq = clean(raw)
    guides = find_guides(seq)
    for g, n in zip(guides, count_offtargets(seq, guides)):
        g["gc"] = gc_content(g["guide"])
        g["offtargets"] = n
        g["flags"] = [f for f, bad in (("poly-T", "TTTT" in g["guide"]),
                                       ("GC outside 40-70%", not 0.4 <= g["gc"] <= 0.7)) if bad]
        g["risk"] = "low" if n == 0 and not g["flags"] else ("high" if n > 2 else "medium")
    return seq, guides

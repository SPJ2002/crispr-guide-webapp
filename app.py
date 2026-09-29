import sqlite3
import time
from pathlib import Path
from flask import Flask, jsonify, render_template, request
from core import analyze

MAX_LEN = 10000
app = Flask(__name__)
app.config["DB"] = str(Path(__file__).with_name("guides.db"))


def db():
    conn = sqlite3.connect(app.config["DB"])
    conn.execute("""CREATE TABLE IF NOT EXISTS runs (id INTEGER PRIMARY KEY, created TEXT DEFAULT CURRENT_TIMESTAMP,
                    seq_len INTEGER, n_guides INTEGER, n_low_risk INTEGER, seconds REAL)""")
    return conn


@app.get("/")
def index():
    return render_template("index.html", max_len=MAX_LEN)


@app.post("/api/analyze")
def api_analyze():
    raw = (request.get_json(silent=True) or {}).get("sequence", "")
    t0 = time.perf_counter()
    try:
        seq, guides = analyze(raw)
    except ValueError as e:
        return jsonify(error=str(e)), 400
    if len(seq) > MAX_LEN:
        return jsonify(error=f"Sequence too long (max {MAX_LEN} bp)"), 413
    secs = round(time.perf_counter() - t0, 3)
    low = sum(g["risk"] == "low" for g in guides)
    with db() as conn:
        conn.execute("INSERT INTO runs (seq_len, n_guides, n_low_risk, seconds) VALUES (?,?,?,?)",
                     (len(seq), len(guides), low, secs))
    return jsonify(seq_len=len(seq), n_guides=len(guides), n_low_risk=low, seconds=secs, guides=guides)


@app.get("/api/history")
def api_history():
    with db() as conn:
        rows = conn.execute("SELECT id, created, seq_len, n_guides, n_low_risk, seconds FROM runs "
                            "ORDER BY id DESC LIMIT 10").fetchall()
    keys = ["id", "created", "seq_len", "n_guides", "n_low_risk", "seconds"]
    return jsonify([dict(zip(keys, r)) for r in rows])


if __name__ == "__main__":
    app.run(debug=True)

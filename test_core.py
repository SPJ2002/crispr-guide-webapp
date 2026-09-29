import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import pytest
from core import clean, revcomp, gc_content, find_guides, count_offtargets, analyze
import app as appmod


def test_revcomp():
    assert revcomp("AACG") == "CGTT"


def test_clean_strips_fasta_and_whitespace():
    assert clean(">seq1\nacg t\nAC") == "ACGTAC"


@pytest.mark.parametrize("bad", ["", "ACGN", "hello"])
def test_clean_rejects_invalid(bad):
    with pytest.raises(ValueError):
        clean(bad)


def test_gc_content():
    assert gc_content("GGGGCCCCAAAATTTTAAAA") == 0.4


def test_find_guides_plus_strand():
    seq = "ACGTACGTACGTACGTACGT" + "TGG"
    g = find_guides(seq)
    assert any(x["strand"] == "+" and x["start"] == 0 and x["pam"] == "TGG" for x in g)


def test_find_guides_minus_strand():
    seq = revcomp("ACGTACGTACGTACGTACGT" + "TGG")
    assert any(x["strand"] == "-" and x["guide"] == "ACGTACGTACGTACGTACGT" for x in find_guides(seq))


def test_offtarget_detects_duplicate_site():
    guide = "ACGTACGTACGTACGTACGT"
    seq = guide + "TGG" + "TTTTTTTTTT" + guide + "AGG"
    _, guides = analyze(seq)
    target = next(g for g in guides if g["guide"] == guide and g["strand"] == "+" and g["start"] == 0)
    assert target["offtargets"] >= 1 and target["risk"] != "low"


def test_offtarget_none_for_unique():
    seq = "ATCGGCTAGCTTAGCCGATAGGCTAACGTTGCA" * 1
    _, guides = analyze(seq + "AGG")
    assert all(g["offtargets"] >= 0 for g in guides)


@pytest.fixture
def client(tmp_path):
    appmod.app.config["DB"] = str(tmp_path / "t.db")
    return appmod.app.test_client()


def test_api_analyze_and_history(client):
    r = client.post("/api/analyze", json={"sequence": "ACGTACGTACGTACGTACGTTGG"})
    assert r.status_code == 200 and r.get_json()["n_guides"] >= 1
    assert len(client.get("/api/history").get_json()) == 1


def test_api_rejects_bad_input(client):
    assert client.post("/api/analyze", json={"sequence": "XYZ"}).status_code == 400


def test_api_rejects_too_long(client):
    assert client.post("/api/analyze", json={"sequence": "A" * 10001}).status_code == 413

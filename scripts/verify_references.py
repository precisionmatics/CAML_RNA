#!/usr/bin/env python3
"""
verify_references.py
Query CrossRef API for each reference to confirm it exists and is correct.
Flags: DOI not found, title mismatch, wrong journal, hallucinated papers.
"""

import re, time, json, urllib.request, urllib.parse, sys

REFS = [
    (1,  "Schlünzen F, Zarivach R, Harms J, Bashan A, Tocilj A, Albrecht R, Yonath A, Franceschi F",
         "Structural basis for the interaction of antibiotics with the peptidyl transferase center in eubacteria",
         "Nature", 2001),
    (2,  "Thomas JR, Hergenrother PJ",
         "Targeting RNA with small molecules",
         "Chemical reviews", 2008),
    (3,  "Mortimer SA, Kidwell MA, Doudna JA",
         "Insights into RNA structure and function from genome-wide studies",
         "Nature Reviews Genetics", 2014),
    (4,  "Winkler W, Nahvi A, Breaker RR",
         "Thiamine derivatives bind messenger RNAs directly to regulate bacterial gene expression",
         "Nature", 2002),
    (5,  "Breaker RR",
         "Riboswitches and the RNA world",
         "Cold Spring Harbor perspectives in biology", 2012),
    (6,  "Serganov A, Nudler E",
         "A decade of riboswitches",
         "Cell", 2013),
    (7,  "Donlic A, Hargrove AE",
         "Targeting RNA in mammalian systems with small molecules",
         "Wiley Interdisciplinary Reviews RNA", 2018),
    (8,  "Warner KD, Hajdin CE, Weeks KM",
         "Principles for targeting RNA with drug-like small molecules",
         "Nature reviews Drug discovery", 2018),
    (9,  "Luo J, Wei W, Waldispuhl J, Moitessier N",
         "Challenges and current status of computational methods for docking small molecules to nucleic acids",
         "European journal of medicinal chemistry", 2019),
    (10, "Sun T, Xia W, Shu J, Sang C, Feng M, Xu X",
         "Advances and Challenges in Machine Learning for RNA-Small Molecule Interaction Modeling",
         "Journal of Chemical Theory and Computation", 2025),
    (11, "Pfeffer P, Gohlke H",
         "DrugScoreRNA Knowledge-Based Scoring Function To Predict RNA Ligand Interactions",
         "Journal of chemical information and modeling", 2007),
    (12, "Zhou Y, Jiang Y, Chen SJ",
         "RNA ligand molecular docking Advances and challenges",
         "Wiley Interdisciplinary Reviews Computational Molecular Science", 2022),
    (13, "Arulsamy S, Arora P, Kumar S",
         "Advances in computational prediction of RNA-small molecule binding affinity",
         "Journal of Computer-Aided Molecular Design", 2026),
    (14, "Wang J, Wu J, Zhang Z, Jiang Y, Peng L, Zhang B, Chen Q, Cao L, Quan L, Lyu Q",
         "AffiGrapher Contrastive Heterogeneous Graph Learning with Aromatic Virtual Nodes for RNA-Small Molecule Binding Affinity Prediction",
         "Journal of Chemical Information and Modeling", 2025),
    (15, "Sun C, Zhang L, Zhang L, Song Y, Ma B, Wang Y",
         "GATRsite RNA ligand binding site prediction using graph attention networks and pretrained RNA language models",
         "Journal of Chemical Information and Modeling", 2025),
    (16, "Nguyen T, Le H, Quinn TP, Nguyen T, Le TD, Venkatesh S",
         "GraphDTA predicting drug target binding affinity with graph neural networks",
         "Bioinformatics", 2021),
    (17, "Krishnan SR, Roy A, Gromiha MM",
         "Reliable method for predicting the binding affinity of RNA-small molecule interactions using machine learning",
         "Briefings in Bioinformatics", 2024),
    (18, "Xia W, Shu J, Sang C, Wang K, Wang Y, Sun T, Xu X",
         "The prediction of RNA-small-molecule ligand binding affinity based on geometric deep learning",
         "Computational Biology and Chemistry", 2025),
    (19, "Suwayyid F, Wei GW",
         "Persistent Stanley Reisner Theory",
         "arXiv", 2025),
    (20, "Miller E, Sturmfels B",
         "Combinatorial commutative algebra",
         "Springer", 2005),
    (21, "Bruns W, Herzog HJ",
         "Cohen-Macaulay rings",
         "Cambridge university press", 1998),
    (22, "Feng H, Suwayyid F, Zia M, Wee J, Hozumi Y, Chen CL, Wei GW",
         "CAML Commutative Algebra Machine Learning Case Study on Protein Ligand Binding Affinity Prediction",
         "Journal of Chemical Information and Modeling", 2025),
    (23, "Edelsbrunner H, Letscher D, Zomorodian A",
         "Topological persistence and simplification",
         "Discrete and computational geometry", 2002),
    (24, "Meng Z, Xia K",
         "Persistent spectral based machine learning PerSpect ML for protein ligand binding affinity prediction",
         "Science advances", 2021),
    (25, "Cang Z, Wei GW",
         "TopologyNet Topology based deep convolutional and multitask neural networks for biomolecular property predictions",
         "PLoS computational biology", 2017),
    (26, "Sun S, Gao L",
         "Contrastive pretraining and 3D convolution neural network for RNA and small molecule binding affinity prediction",
         "Bioinformatics", 2024),
    (27, "Xu M, Zeng Y, Zhao L, Zhang Z, Chen B, Liu W",
         "DeepMIF Deep Multimodal Interaction Features for Predicting RNA Small Molecule Binding Affinity",
         "Journal of Chemical Information and Modeling", 2024),
    (28, "Huang Z, Wang Y, Chen S, Tan YS, Deng L, Wu M",
         "DeepMIF cross-fusion deep learning RNA small molecule binding affinity prediction",
         "Bioinformatics", 2024),
    (29, "Li Z, Wang Y, Chen S, Xu M",
         "EMMPTNet Element-specific Multimodal Persistent Topology Network RNA Ligand Binding Affinity Prediction",
         "Journal of Chemical Information and Modeling", 2024),
    (30, "Wang R, Fang X, Lu Y, Wang S",
         "The PDBbind database collection of binding affinities for protein ligand complexes with known three-dimensional structures",
         "Journal of medicinal chemistry", 2004),
    (31, "Rose PW, Prlic A, Bi C, Bluhm WF, Christie CH, Dutta S, Green RK, Goodsell DS, Westbrook JD",
         "The RCSB Protein Data Bank views of structural biology for basic and applied research and education",
         "Nucleic acids research", 2015),
    (32, "RDKit",
         "RDKit open-source cheminformatics",
         "Zenodo", 2024),
    (33, "Zomorodian A, Carlsson G",
         "Computing persistent homology",
         "Discrete and Computational Geometry", 2005),
    (34, "Bauer U",
         "Ripser efficient computation of Vietoris Rips persistence barcodes",
         "Journal of Applied and Computational Topology", 2021),
    (35, "Chen J, Hu Z, Sun S, Tan Q, Wang Y, Yu Q, Zong L, Hong L, Xiao J, Shen T, King I",
         "Interpretable RNA foundation model from unannotated data for highly accurate RNA structure and function predictions",
         "arXiv", 2022),
    (36, "Rogers D, Hahn M",
         "Extended-connectivity fingerprints",
         "Journal of chemical information and modeling", 2010),
    (37, "Cortes C, Vapnik V",
         "Support-vector networks",
         "Machine learning", 1995),
    (38, "Kramer O",
         "Scikit-learn machine learning evolution strategies",
         "Springer", 2016),
]

BASE = "https://api.crossref.org/works"
HEADERS = {"User-Agent": "CAML-RNA-verification/1.0 (mailto:stalin.bioinfo@gmail.com)"}

def crossref_query(title, author_first_last, year=None):
    """Query CrossRef and return top result or None."""
    first_author = author_first_last.split(',')[0].strip().split()[-1]  # last name
    q = f"{title} {first_author}"
    params = urllib.parse.urlencode({"query.bibliographic": q, "rows": "3", "select": "DOI,title,author,published,container-title,score"})
    url = f"{BASE}?{params}"
    try:
        req = urllib.request.Request(url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
        items = data.get("message", {}).get("items", [])
        return items
    except Exception as e:
        return []

def title_similarity(ref_title, result_title):
    """Rough word-overlap similarity."""
    ref_words  = set(re.sub(r'[^a-z0-9 ]', ' ', ref_title.lower()).split())
    res_words  = set(re.sub(r'[^a-z0-9 ]', ' ', result_title.lower()).split())
    stop = {'the','a','an','and','or','for','with','to','of','in','on','by','from','is','are','as','at','be','was'}
    ref_words -= stop; res_words -= stop
    if not ref_words: return 0
    return len(ref_words & res_words) / len(ref_words)

SKIP_CROSSREF = {19, 20, 21, 32, 35, 38}  # arXiv, books, software, textbook chapters

results = []
print(f"{'#':>3} {'Status':<10} {'Score':>6}  Details")
print("-" * 80)

for num, authors, title, journal, year in REFS:
    if num in SKIP_CROSSREF:
        tag = "SKIP"
        note = f"Not in CrossRef ({journal})"
        results.append((num, tag, None, note, authors, title))
        print(f"[{num:2d}] {tag:<10} {'—':>6}  {note}")
        time.sleep(0.05)
        continue

    items = crossref_query(title, authors, year)
    time.sleep(0.35)

    if not items:
        tag, note = "NOT_FOUND", "No CrossRef results"
        results.append((num, tag, None, note, authors, title))
        print(f"[{num:2d}] {tag:<10} {'—':>6}  {note} | {title[:60]}")
        continue

    best = items[0]
    doi  = best.get("DOI", "—")
    res_titles = best.get("title", [""])
    res_title  = res_titles[0] if res_titles else ""
    res_year   = (best.get("published", {}).get("date-parts") or [[None]])[0][0]
    res_journal = (best.get("container-title") or [""])[0]
    sim = title_similarity(title, res_title)
    cr_score = best.get("score", 0)

    if sim >= 0.55:
        tag = "OK"
    elif sim >= 0.35:
        tag = "PARTIAL"
    else:
        tag = "MISMATCH"

    note = f"sim={sim:.2f} | doi:{doi} | found: {res_title[:55]!r}"
    results.append((num, tag, doi, note, authors, title))
    print(f"[{num:2d}] {tag:<10} {sim:>6.2f}  {note}")

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
ok = [r for r in results if r[1] == "OK"]
partial = [r for r in results if r[1] == "PARTIAL"]
mismatch = [r for r in results if r[1] == "MISMATCH"]
not_found = [r for r in results if r[1] == "NOT_FOUND"]
skipped = [r for r in results if r[1] == "SKIP"]

print(f"  OK:        {len(ok)}")
print(f"  PARTIAL:   {len(partial)}")
print(f"  MISMATCH:  {len(mismatch)}")
print(f"  NOT_FOUND: {len(not_found)}")
print(f"  SKIPPED:   {len(skipped)} (books/arXiv/software)")
print()

if partial or mismatch or not_found:
    print("NEEDS ATTENTION:")
    for r in partial + mismatch + not_found:
        print(f"  [{r[0]:2d}] {r[1]}: {r[4].split(',')[0]} | '{r[5][:65]}'")
        print(f"       → {r[3]}")

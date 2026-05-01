"""
nr_cofactors.py
---------------
Structural bioinformatics pipeline for analysing nuclear receptor (NR)
cofactor complexes from the Protein Data Bank (PDB).

Research conducted at Nazarbayev University, BIOL363, under Dr. Ferdinand Molnár.

Pipeline overview:
    1. Fetch UniProt IDs for 48 human nuclear receptor genes via UniProt REST API
    2. Retrieve protein details and extract ligand-binding domain (LBD) sequences
    3. Fetch associated PDB structure IDs for each NR via RCSB Search API
    4. Download PDB XML files and classify structures as DBD (contains Zinc) or LBD
    5. Filter LBD structures for NR–cofactor peptide complexes
    6. Visualise results

Usage:
    python nr_cofactors.py

All outputs are written to OUTPUT_DIR.
"""

import os
import re
import sys
import shutil
import textwrap
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from time import sleep
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import requests

# ─────────────────────────────────────────────
# CONFIGURATION — edit these paths before running
# ─────────────────────────────────────────────

BASE_DIR   = Path("/path/to/biol363_project")
OUTPUT_DIR = BASE_DIR / "outputs"

# Intermediate data directories
ALL_NR_STRUCTURES_DIR = OUTPUT_DIR / "all_nr_structures"
LBD_STRUCTURES_DIR    = OUTPUT_DIR / "lbd_structures"
PIPELINE_DATA_DIR     = OUTPUT_DIR / "nr_pipeline_data"

# API endpoints
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL  = "https://rest.uniprot.org/uniprotkb/{uid}"
PDB_SEARCH_URL     = "https://search.rcsb.org/rcsbsearch/v2/query"
PDB_DOWNLOAD_URL   = "https://files.rcsb.org/download/{pdb_id}.xml"

# Filtering thresholds
MIN_RECEPTOR_LENGTH = 50   # amino acids
MAX_PEPTIDE_LENGTH  = 49   # amino acids

# API rate limiting
API_SLEEP = 0.5  # seconds between requests

# Known cofactor UniProt IDs
KNOWN_COFACTORS = {
    "Q15788": "SRC-1", "Q15596": "SRC-2", "Q9Y6Q9": "SRC-3",
    "Q92793": "CBP",   "Q09472": "p300",  "O75376": "NCoR",  "Q9Y618": "SMRT",
}

# 48 official HGNC gene symbols for human nuclear receptors
HUMAN_NUCLEAR_RECEPTOR_GENES = sorted([
    "AR",    "ESR1",  "ESR2",  "ESRRA", "ESRRB", "ESRRG", "HNF4A", "HNF4G",
    "NR0B1", "NR0B2", "NR1D1", "NR1D2", "NR1H2", "NR1H3", "NR1H4", "NR1I2",
    "NR1I3", "VDR",   "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2",
    "NR2F6", "NR3C1", "NR3C2", "PGR",   "NR4A1", "NR4A2", "NR4A3", "NR5A1",
    "NR5A2", "NR6A1", "PPARA", "PPARD", "PPARG", "RARA",  "RARB",  "RARG",
    "RORA",  "RORB",  "RORC",  "RXRA",  "RXRB",  "RXRG",  "THRA",  "THRB",
])


# ─────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────

def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def safe_text(elem, default: str = "") -> str:
    """Return stripped text from an XML element, or default if absent."""
    return elem.text.strip() if elem is not None and elem.text else default


# ─────────────────────────────────────────────
# STEP 1 — FETCH UNIPROT IDs
# ─────────────────────────────────────────────

def get_uniprot_id(gene_symbol: str) -> Optional[str]:
    """
    Query UniProt for the reviewed human entry matching a gene symbol.
    Returns the primary accession number, or None if not found.
    """
    params = {
        "query":  f"(gene:{gene_symbol}) AND (organism_id:9606) AND (reviewed:true)",
        "format": "json",
        "fields": "accession",
        "size":   1,
    }
    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=15)
        response.raise_for_status()
        results = response.json().get("results", [])
        return results[0]["primaryAccession"] if results else None
    except requests.exceptions.RequestException as e:
        print(f"  Error fetching UniProt ID for {gene_symbol}: {e}", file=sys.stderr)
        return None


def fetch_all_uniprot_ids(genes: List[str]) -> Tuple[dict, List[str]]:
    """
    Fetch UniProt IDs for a list of gene symbols.
    Returns (gene → uid mapping, list of genes with no match).
    """
    print("Step 1: Fetching UniProt IDs for human nuclear receptors...\n")
    id_map, not_found = {}, []

    for i, gene in enumerate(genes):
        print(f"  [{i+1}/{len(genes)}] {gene}...", end=" ")
        uid = get_uniprot_id(gene)
        if uid:
            id_map[gene] = uid
            print(uid)
        else:
            not_found.append(gene)
            print("not found")
        sleep(API_SLEEP)

    print(f"\nFound {len(id_map)}/{len(genes)} UniProt IDs.")
    return id_map, not_found


# ─────────────────────────────────────────────
# STEP 2 — RETRIEVE PROTEIN DETAILS & EXTRACT LBDs
# ─────────────────────────────────────────────

def get_protein_details_and_lbd(uniprot_id: str) -> Optional[dict]:
    """
    Fetch protein details from UniProt for a given accession.
    Extracts protein name, full sequence length, and annotated LBD coordinates.
    Also checks for presence of zinc finger domains (indicative of DBD).
    """
    url    = UNIPROT_ENTRY_URL.format(uid=uniprot_id)
    params = {"format": "json", "fields": "protein_name,sequence,ft_domain,ft_zn_fing"}
    try:
        response = requests.get(url, params=params, timeout=15)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"  Error fetching details for {uniprot_id}: {e}", file=sys.stderr)
        return None

    name     = (data.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", "N/A"))
    sequence = data.get("sequence", {}).get("value", "")
    length   = data.get("sequence", {}).get("length", 0)

    lbd = {"found": False, "location": "Not found", "sequence": ""}
    has_zn_finger = False

    for feature in data.get("features", []):
        desc = feature.get("description", {})
        desc_text = desc.get("value", "") if isinstance(desc, dict) else str(desc)

        if feature.get("type") == "Domain" and "Ligand-binding" in desc_text and not lbd["found"]:
            loc   = feature.get("location", {})
            start = loc.get("start", {}).get("value")
            end   = loc.get("end", {}).get("value")
            if start and end:
                lbd = {
                    "found":    True,
                    "location": f"{start}-{end}",
                    "sequence": sequence[start - 1:end],
                }

        if feature.get("type") == "Zinc finger":
            has_zn_finger = True

    return {
        "name":           name,
        "length":         length,
        "lbd_found":      lbd["found"],
        "lbd_location":   lbd["location"],
        "lbd_sequence":   lbd["sequence"],
        "has_zn_finger":  has_zn_finger,
    }


def analyze_and_extract_lbds(id_map: dict) -> List[dict]:
    """Retrieve protein details for each unique UniProt ID and extract LBD info."""
    print("\nStep 2: Retrieving protein details and extracting LBD sequences...")
    results     = []
    unique_ids  = sorted(set(id_map.values()))

    for i, uid in enumerate(unique_ids):
        genes    = ", ".join(g for g, u in id_map.items() if u == uid)
        print(f"  [{i+1}/{len(unique_ids)}] {uid} ({genes})...", end=" ")
        details  = get_protein_details_and_lbd(uid)
        sleep(API_SLEEP)

        if details:
            details.update({"gene": genes, "uniprot_id": uid, "status": "OK"})
            print(f"{'LBD found' if details['lbd_found'] else 'no LBD annotated'}")
        else:
            details = {"gene": genes, "uniprot_id": uid, "status": "Details fetch failed"}
            print("failed")
        results.append(details)

    return results


def write_lbd_report(results: List[dict], filename: Path) -> None:
    """Write a tabular summary of LBD extraction results."""
    with open(filename, "w", encoding="utf-8") as f:
        f.write("--- Human Nuclear Receptor LBD Extraction Analysis ---\n\n")
        header = f"{'Gene(s)':<15} | {'UniProt ID':<12} | {'ZN Finger?':<12} | {'LBD Found?':<11} | {'LBD Location':<13} | {'LBD Len':<9} | Protein Name\n"
        f.write(header)
        f.write("-" * (len(header) + 5) + "\n")

        for r in sorted(results, key=lambda x: x["gene"]):
            if r["status"] == "OK":
                f.write(
                    f"{r['gene']:<15} | {r['uniprot_id']:<12} | "
                    f"{'Yes' if r['has_zn_finger'] else 'No':<12} | "
                    f"{'Yes' if r['lbd_found'] else 'No':<11} | "
                    f"{r['lbd_location']:<13} | "
                    f"{len(r['lbd_sequence']):<9} | {r['name']}\n"
                )
            else:
                f.write(f"{r['gene']:<15} | {r['uniprot_id']:<12} | {'N/A':<12} | {'N/A':<11} | {'N/A':<13} | {'N/A':<9} | {r['status']}\n")
    print(f"  LBD report saved → {filename}")


def write_fasta(results: List[dict], filename: Path) -> None:
    """Export extracted LBD sequences to FASTA format."""
    with open(filename, "w", encoding="utf-8") as f:
        for r in sorted(results, key=lambda x: x["gene"]):
            if r["status"] == "OK" and r["lbd_found"]:
                warning = " [WARNING: ZN-FINGER DOMAIN PRESENT]" if r["has_zn_finger"] else ""
                f.write(f">{r['uniprot_id']}|{r['gene']} {r['name']} [LBD: {r['lbd_location']}]{warning}\n")
                f.write(textwrap.fill(r["lbd_sequence"], 80) + "\n")
    print(f"  FASTA file saved → {filename}")


# ─────────────────────────────────────────────
# STEP 3 — FETCH PDB STRUCTURE IDs
# ─────────────────────────────────────────────

def fetch_pdb_ids_for_uniprot(uniprot_id: str) -> Tuple[str, List[str]]:
    """Query RCSB for all PDB entries associated with a UniProt accession."""
    query = {
        "query": {
            "type": "terminal", "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    try:
        response = requests.post(PDB_SEARCH_URL, json=query, timeout=30)
        response.raise_for_status()
        pdb_ids = [item["identifier"] for item in response.json()["result_set"]]
        return uniprot_id, pdb_ids
    except requests.exceptions.RequestException as e:
        print(f"  Error fetching PDB IDs for {uniprot_id}: {e}", file=sys.stderr)
        return uniprot_id, []


def fetch_all_pdb_ids(uniprot_ids: List[str]) -> dict:
    """Fetch PDB IDs for all UniProt IDs in parallel."""
    print("\nStep 3: Fetching PDB IDs for all NRs (parallel)...")
    pdb_map = {}
    with ThreadPoolExecutor(max_workers=5) as executor:
        for uid, pdb_ids in executor.map(fetch_pdb_ids_for_uniprot, uniprot_ids):
            pdb_map[uid] = pdb_ids
            print(f"  {uid}: {len(pdb_ids)} structures")
    return pdb_map


# ─────────────────────────────────────────────
# STEP 4 — DOWNLOAD & CLASSIFY PDB STRUCTURES
# ─────────────────────────────────────────────

def download_pdb_xml(pdb_id: str, directory: Path, max_attempts: int = 3) -> bool:
    """
    Download a PDB XML file into the given directory.
    Skips files that already exist and are valid XML.
    Retries up to max_attempts times with exponential backoff.
    """
    mkdir(directory)
    filepath = directory / f"{pdb_id}.xml"

    if filepath.exists():
        try:
            ET.parse(filepath)
            return True  # already downloaded and valid
        except ET.ParseError:
            print(f"  {pdb_id}.xml is corrupted, re-downloading...")

    for attempt in range(1, max_attempts + 1):
        try:
            url      = PDB_DOWNLOAD_URL.format(pdb_id=pdb_id)
            response = requests.get(url, timeout=20)
            response.raise_for_status()
            filepath.write_bytes(response.content)
            ET.parse(filepath)  # validate
            return True
        except ET.ParseError:
            print(f"  Downloaded {pdb_id}.xml is corrupt (attempt {attempt})")
        except requests.exceptions.RequestException as e:
            print(f"  Failed to download {pdb_id} (attempt {attempt}): {e}")
        if attempt < max_attempts:
            sleep(2 ** attempt)

    return False


def contains_zinc(filepath: Path) -> bool:
    """
    Return True if a PDB XML file contains an atomic Zinc (ZN) entry.
    Zinc presence indicates a DNA-binding domain (DBD) structure.
    """
    ns = {"pdbx": "http://pdbml.pdb.org/schema/pdbx-v50.xsd"}
    try:
        root = ET.parse(filepath).getroot()
        return root.find('.//pdbx:atom_site/pdbx:type_symbol[.="ZN"]', ns) is not None
    except ET.ParseError:
        return False  # treat unparseable files conservatively


def download_and_classify_all(pdb_map: dict) -> List[dict]:
    """
    Download all PDB XML files organised by UniProt ID,
    then classify each as DBD (Zinc present) or LBD (no Zinc).
    Returns a list of per-receptor summary dicts.
    """
    print("\nStep 4: Downloading and classifying PDB structures...")
    download_tasks = []
    for uid, pdb_ids in pdb_map.items():
        nr_dir = ALL_NR_STRUCTURES_DIR / uid
        for pdb_id in pdb_ids:
            download_tasks.append((pdb_id, nr_dir))

    unique_pdbs = {t[0] for t in download_tasks}
    print(f"  Downloading {len(unique_pdbs)} unique structures...")
    with ThreadPoolExecutor(max_workers=10) as executor:
        for i, _ in enumerate(executor.map(lambda t: download_pdb_xml(*t), download_tasks)):
            if i % 100 == 0:
                sys.stdout.write(f"\r  {i+1}/{len(download_tasks)} downloaded...")
                sys.stdout.flush()
    print("\n  Download complete.")

    # Classify each structure
    report_data = []
    for uid, pdb_ids in pdb_map.items():
        nr_dir    = ALL_NR_STRUCTURES_DIR / uid
        dbd_count = sum(contains_zinc(nr_dir / f"{p}.xml") for p in pdb_ids if (nr_dir / f"{p}.xml").exists())
        lbd_count = len(pdb_ids) - dbd_count
        report_data.append({"uniprot_id": uid, "total": len(pdb_ids), "dbd": dbd_count, "lbd": lbd_count})

    return report_data


def write_classification_report(report_data: List[dict], filename: Path) -> None:
    """Write the DBD/LBD classification summary to a text file."""
    with open(filename, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("=== Full Nuclear Receptor Structural Classification Report ===\n")
        f.write("=" * 70 + "\n\n")
        f.write("Classification: structures with atomic Zinc → DBD; without → LBD candidate.\n\n")
        header = f"{'UniProt ID':<12} | {'Total':<8} | {'DBD (w/ ZN)':<14} | {'LBD (no ZN)'}\n"
        f.write(header + "-" * (len(header) + 5) + "\n")
        for item in sorted(report_data, key=lambda x: x["uniprot_id"]):
            f.write(f"{item['uniprot_id']:<12} | {item['total']:<8} | {item['dbd']:<14} | {item['lbd']}\n")
    print(f"  Classification report saved → {filename}")


# ─────────────────────────────────────────────
# STEP 5 — FILTER LBD STRUCTURES & EXTRACT COMPLEXES
# ─────────────────────────────────────────────

def copy_lbd_structures() -> None:
    """Copy all non-Zinc PDB XML files from all_nr_structures into lbd_structures."""
    mkdir(LBD_STRUCTURES_DIR)
    print(f"\nStep 5a: Copying LBD structures → {LBD_STRUCTURES_DIR}")
    total, copied = 0, 0
    for xml_file in ALL_NR_STRUCTURES_DIR.rglob("*.xml"):
        total += 1
        if not contains_zinc(xml_file):
            dest = LBD_STRUCTURES_DIR / xml_file.name
            shutil.copy2(xml_file, dest)
            copied += 1
    print(f"  Scanned {total} files, copied {copied} LBD candidates.")


@dataclass
class Chain:
    entity:  str
    seq:     str
    length:  int
    desc:    str
    uniprot: Optional[str]


def parse_pdb_xml(pdb_id: str) -> Tuple[List[Chain], str]:
    """
    Parse a PDB XML file from LBD_STRUCTURES_DIR.
    Returns a list of polypeptide chains and the structure title.
    """
    filepath = LBD_STRUCTURES_DIR / f"{pdb_id}.xml"
    if not filepath.exists():
        return [], f"{pdb_id}: XML missing"

    ns = {"pdbx": "http://pdbml.pdb.org/schema/pdbx-v50.xsd"}
    try:
        root   = ET.parse(filepath).getroot()
        title  = safe_text(root.find(".//pdbx:struct_title", ns))
        chains = []

        for poly in root.findall(".//pdbx:entity_poly", ns):
            if safe_text(poly.find("pdbx:type", ns)) != "polypeptide(L)":
                continue
            eid = safe_text(poly.find("pdbx:entity_id", ns))
            if not eid:
                continue

            seq_el = (poly.find("pdbx:pdbx_seq_one_letter_code_can", ns)
                      or poly.find("pdbx:pdbx_seq_one_letter_code", ns))
            seq = safe_text(seq_el).replace("\n", "").replace(" ", "")
            if not seq:
                continue

            desc   = safe_text(root.find(f".//pdbx:entity[@id='{eid}']/pdbx:pdbx_description", ns), "N/A")
            up_el  = root.find(f".//pdbx:entity[@id='{eid}']/pdbx:db_reference[@db='UniProt']", ns)
            uniprot = up_el.get("id") if up_el is not None else None
            chains.append(Chain(eid, seq, len(seq), desc, uniprot))

        return chains, title
    except ET.ParseError as e:
        return [], f"{pdb_id}: parse error – {e}"


def process_pdb_for_complex(pdb_id: str) -> Optional[str]:
    """
    Determine whether a PDB structure qualifies as an NR–cofactor peptide complex:
    - At least one receptor chain (≥ MIN_RECEPTOR_LENGTH AA)
    - At least one peptide chain (≤ MAX_PEPTIDE_LENGTH AA) that is not
      a subsequence of any receptor chain

    Returns a formatted detail string if it qualifies, otherwise None.
    """
    chains, _ = parse_pdb_xml(pdb_id)
    if not chains:
        return None

    receptors = [c for c in chains if c.length >= MIN_RECEPTOR_LENGTH]
    peptides  = [c for c in chains if c.length <= MAX_PEPTIDE_LENGTH]

    if not receptors or not peptides:
        return None

    valid_peptides = [p for p in peptides if not any(p.seq in r.seq for r in receptors)]
    if not valid_peptides:
        return None

    lines = [f"{pdb_id}:"]
    for r in receptors:
        cof = KNOWN_COFACTORS.get(r.uniprot, "Receptor")
        lines.append(f"  - Entity {r.entity}: {r.length} AA | {r.desc} | {r.uniprot or '-'} ({cof})")
    for p in valid_peptides:
        cof = KNOWN_COFACTORS.get(p.uniprot, "Unknown")
        lines.append(f"  - Entity {p.entity}: {p.length} AA | {p.desc} | {p.uniprot or '-'} ({cof})")
    return "\n".join(lines)


def filter_nr_peptide_complexes() -> List[str]:
    """Run parallel filtering of LBD structures for NR–cofactor complexes."""
    mkdir(PIPELINE_DATA_DIR)
    xml_files = [f.stem.upper() for f in LBD_STRUCTURES_DIR.glob("*.xml")]
    print(f"\nStep 5b: Filtering {len(xml_files)} LBD structures for NR–peptide complexes...")

    results = []
    with ThreadPoolExecutor(max_workers=16) as executor:
        for i, res in enumerate(executor.map(process_pdb_for_complex, xml_files)):
            if res:
                results.append(res)
            if i % 100 == 0:
                print(f"\r  Processed {i+1}/{len(xml_files)} | Kept {len(results)}", end="")
    print()

    complex_ids = [r.splitlines()[0][:-1] for r in results]
    (PIPELINE_DATA_DIR / "nr_peptide_complexes.txt").write_text("\n".join(complex_ids))
    (PIPELINE_DATA_DIR / "nr_peptide_complexes_detailed.txt").write_text("\n\n".join(results))
    print(f"  Found {len(complex_ids)} NR–peptide complexes.")
    print(f"  Results saved → {PIPELINE_DATA_DIR}")
    return complex_ids


# ─────────────────────────────────────────────
# STEP 6 — VISUALISATION
# ─────────────────────────────────────────────

def plot_dbd_lbd_distribution(report_data: List[dict], output_path: Path) -> None:
    """Stacked bar chart of DBD vs LBD structure counts per nuclear receptor."""
    df = pd.DataFrame(report_data).set_index("uniprot_id")[["dbd", "lbd"]]
    df.columns = ["DBD (w/ ZN)", "LBD (no ZN)"]

    fig, ax = plt.subplots(figsize=(16, 6))
    df.plot(kind="bar", stacked=True, ax=ax, color={"DBD (w/ ZN)": "orange", "LBD (no ZN)": "steelblue"})
    ax.set_title("DBD vs LBD Structures per Nuclear Receptor", fontsize=14)
    ax.set_xlabel("UniProt ID")
    ax.set_ylabel("Number of PDB Structures")
    ax.legend(title="Domain")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.show()
    print(f"  Plot saved → {output_path}")


def plot_chain_distribution(detailed_file: Path, output_path: Path) -> None:
    """Bar chart showing how many unique polypeptide chains are in each complex."""
    pdb_counts = {}
    current_pdb = None

    for line in detailed_file.read_text().splitlines():
        line = line.strip()
        if line.endswith(":") and "#" not in line:
            current_pdb = line[:-1]
            pdb_counts[current_pdb] = 0
        elif line.startswith("- Entity") and current_pdb:
            pdb_counts[current_pdb] += 1

    multi_chain = sum(1 for c in pdb_counts.values() if c > 2)
    print(f"  Complexes with >2 polypeptide chains: {multi_chain}")

    summary = pd.Series(list(pdb_counts.values())).value_counts().sort_index()
    summary = summary[summary.index >= 2]

    fig, ax = plt.subplots(figsize=(8, 5))
    summary.plot(kind="bar", ax=ax, color="seagreen")
    ax.set_title("Polypeptide Chain Count Distribution in NR–Peptide Complexes")
    ax.set_xlabel("Number of Unique Polypeptide Chains")
    ax.set_ylabel("Number of Complexes")
    for i, v in enumerate(summary):
        ax.text(i, v + 0.1, str(v), ha="center")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.show()
    print(f"  Plot saved → {output_path}")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main() -> None:
    mkdir(OUTPUT_DIR)

    # Step 1: UniProt IDs
    id_map, not_found = fetch_all_uniprot_ids(HUMAN_NUCLEAR_RECEPTOR_GENES)

    # Step 2: LBD sequences
    lbd_results = analyze_and_extract_lbds(id_map)
    write_lbd_report(lbd_results, OUTPUT_DIR / "nuclear_receptor_LBD_analysis.txt")
    write_fasta(lbd_results, OUTPUT_DIR / "nuclear_receptor_LBDs.fasta")

    # Step 3: PDB IDs
    pdb_map = fetch_all_pdb_ids(list(set(id_map.values())))

    # Step 4: Download + classify
    report_data = download_and_classify_all(pdb_map)
    write_classification_report(report_data, OUTPUT_DIR / "nr_structural_classification.txt")

    # Step 5: Filter complexes
    copy_lbd_structures()
    filter_nr_peptide_complexes()

    # Step 6: Visualise
    print("\nStep 6: Generating visualisations...")
    plot_dbd_lbd_distribution(report_data, OUTPUT_DIR / "dbd_vs_lbd.png")
    detailed_file = PIPELINE_DATA_DIR / "nr_peptide_complexes_detailed.txt"
    if detailed_file.exists():
        plot_chain_distribution(detailed_file, OUTPUT_DIR / "chain_distribution.png")

    print("\n=== Pipeline complete ===")
    print(f"All outputs in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

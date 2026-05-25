"""
UniProt API helpers and FASTA utilities.

Provides functions to query protein sequences and metadata from UniProt,
and to write / validate FASTA files in UniProt format.
"""

import logging
import os
import sys
import time
from typing import List, Optional

import requests

logger = logging.getLogger(__name__)

WEBSITE_API = "https://rest.uniprot.org/"
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

_FIELDS = "accession,sequence,protein_name,gene_names,organism_name,organism_id,reviewed,id"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_url(url: str, **kwargs) -> requests.Response:
    response = requests.get(url, **kwargs)
    if not response.ok:
        logger.error("Request failed: %s", response.text)
        response.raise_for_status()
        sys.exit()
    return response


def _build_fasta_header(entry: dict) -> str:
    """Return a UniProt-style FASTA header string for *entry*."""
    prot_name = "Unknown protein"
    desc = entry.get("proteinDescription", {})
    if "recommendedName" in desc and desc["recommendedName"]:
        fn = desc["recommendedName"].get("fullName", "")
        prot_name = fn.get("value", prot_name) if isinstance(fn, dict) else fn or prot_name
    elif "submittedName" in desc and desc["submittedName"]:
        submitted = desc["submittedName"]
        if isinstance(submitted, list) and submitted:
            fn = submitted[0].get("fullName", "")
            prot_name = fn.get("value", prot_name) if isinstance(fn, dict) else fn or prot_name

    accession = entry.get("primaryAccession", "UNKNOWN")
    entry_name = str(entry.get("uniProtkbId") or entry.get("uniProtKBId") or f"{accession}_HUMAN")

    reviewed = entry.get("reviewed", None)
    is_reviewed = reviewed if isinstance(reviewed, bool) else "Swiss-Prot" in entry.get("entryType", "")
    db_type = "sp" if is_reviewed else "tr"

    gene_name = ""
    genes = entry.get("genes", [])
    if isinstance(genes, list) and genes:
        gi = genes[0]
        if isinstance(gi, dict):
            gene_name = gi.get("geneName", {}).get("value", "") or gi.get("geneName", "")

    taxid = entry.get("organism", {}).get("taxonId") if isinstance(entry.get("organism"), dict) else None
    pe = "1" if is_reviewed else "2"
    sv = str(entry.get("version", 1))

    parts = [
        f"{db_type}|{accession}|{entry_name}",
        prot_name.replace("\n", " ").replace("\r", " ").strip(),
        "OS=Homo sapiens",
        f"OX={taxid or 9606}",
    ]
    if gene_name:
        parts.append(f"GN={gene_name}")
    parts += [f"PE={pe}", f"SV={sv}"]
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Batch sequence retrieval
# ---------------------------------------------------------------------------


def get_protein_sequences_batch(
    protein_ids: List[str],
    batch_size: int = 50,
) -> List[dict]:
    """
    Retrieve protein sequences from UniProt in batches.

    Args:
        protein_ids: List of UniProt accession strings.
        batch_size: Number of accessions per request. Default ``50``.

    Returns:
        List of UniProt entry dicts (JSON ``results`` items).
    """
    all_sequences: List[dict] = []
    clean_ids = [pid.strip() for pid in protein_ids if pid and pid.strip()]
    n_batches = (len(clean_ids) - 1) // batch_size + 1

    logger.info("Querying %d proteins from UniProt", len(clean_ids))

    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        logger.info("Batch %d/%d (%d proteins)", batch_num, n_batches, len(batch))

        accessions = " OR ".join(f"accession:{acc}" for acc in batch)
        url = f"{WEBSITE_API}uniprotkb/search?query=({accessions})&fields={_FIELDS}&format=json"

        try:
            data = _get_url(url).json()
            results = data.get("results", [])
            all_sequences.extend(results)
            logger.info("  Retrieved %d sequences", len(results))
        except Exception as exc:
            logger.warning("  Batch %d failed: %s", batch_num, exc)
            continue

        time.sleep(0.1)

    logger.info("Total sequences retrieved: %d", len(all_sequences))
    return all_sequences


def get_swissprot_sequences_batch(
    protein_ids: List[str],
    batch_size: int = 50,
) -> List[dict]:
    """
    Retrieve only Swiss-Prot reviewed sequences from UniProt in batches.

    Args:
        protein_ids: List of UniProt accession strings.
        batch_size: Number of accessions per request. Default ``50``.

    Returns:
        List of Swiss-Prot entry dicts.
    """
    all_sequences: List[dict] = []
    clean_ids = [pid.strip() for pid in protein_ids if pid and pid.strip()]
    n_batches = (len(clean_ids) - 1) // batch_size + 1

    logger.info("Querying %d proteins (Swiss-Prot only)", len(clean_ids))

    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        logger.info("Batch %d/%d (%d proteins)", batch_num, n_batches, len(batch))

        accessions = " OR ".join(f"accession:{acc}" for acc in batch)
        query = f"({accessions}) AND reviewed:true"
        url = f"{WEBSITE_API}uniprotkb/search?query={query}&fields={_FIELDS}&format=json"

        try:
            data = _get_url(url).json()
            reviewed = [e for e in data.get("results", []) if "Swiss-Prot" in e.get("entryType", "")]
            all_sequences.extend(reviewed)
            logger.info("  Retrieved %d Swiss-Prot sequences", len(reviewed))
        except Exception as exc:
            logger.warning("  Batch %d failed: %s", batch_num, exc)
            continue

        time.sleep(0.1)

    logger.info("Total Swiss-Prot sequences retrieved: %d", len(all_sequences))
    return all_sequences


# ---------------------------------------------------------------------------
# Full human proteome
# ---------------------------------------------------------------------------


def query_human_proteome(reviewed_only: bool = True) -> List[dict]:
    """
    Stream the full human proteome from UniProt.

    Args:
        reviewed_only: If ``True``, return only Swiss-Prot reviewed proteins (default).

    Returns:
        List of UniProt entry dicts.
    """
    if reviewed_only:
        query = "organism_id:9606 AND reviewed:true"
    else:
        query = "organism_id:9606"

    count_url = f"{WEBSITE_API}uniprotkb/search?query={query}&size=0"
    try:
        total = int(_get_url(count_url).headers.get("x-total-results", 0))
        logger.info("Human proteins available: %d", total)
    except Exception as exc:
        logger.error("Error fetching count: %s", exc)
        return []

    if total == 0:
        logger.warning("No proteins found for query: %s", query)
        return []

    stream_url = f"{WEBSITE_API}uniprotkb/stream?query={query}&fields={_FIELDS}&format=json"
    logger.info("Downloading %d proteins...", total)

    try:
        entries = _get_url(stream_url).json().get("results", [])
        logger.info("Retrieved %d sequences", len(entries))
        return entries
    except Exception as exc:
        logger.error("Error streaming human proteome: %s", exc)
        return []


# ---------------------------------------------------------------------------
# FASTA export
# ---------------------------------------------------------------------------


def write_fasta(entries: List[dict], path: str) -> int:
    """
    Write a list of UniProt entry dicts to a FASTA file.

    Args:
        entries: UniProt entry dicts (from :func:`get_protein_sequences_batch` etc.).
        path: Output file path.

    Returns:
        Number of proteins written.
    """
    written = 0
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for entry in entries:
            sequence = entry.get("sequence", {}).get("value", "")
            if not sequence:
                continue
            fh.write(f">{_build_fasta_header(entry)}\n")
            for i in range(0, len(sequence), 60):
                fh.write(sequence[i:i + 60] + "\n")
            written += 1
    logger.info("FASTA written: %s (%d proteins)", path, written)
    return written


def create_human_proteome_fasta(
    reviewed_only: bool = True,
    file_loc: str = "",
) -> Optional[str]:
    """
    Download the human proteome from UniProt and write it to a FASTA file.

    Args:
        reviewed_only: If ``True``, include only Swiss-Prot reviewed proteins (default).
        file_loc: Output directory (default: current directory).

    Returns:
        Path to the created FASTA file, or ``None`` on failure.
    """
    file_loc = os.path.normpath(file_loc.strip()) if file_loc.strip() else ""
    if file_loc:
        try:
            os.makedirs(file_loc, exist_ok=True)
        except Exception as exc:
            logger.warning("Could not create '%s': %s. Writing to current directory.", file_loc, exc)
            file_loc = ""

    entries = query_human_proteome(reviewed_only=reviewed_only)
    if not entries:
        logger.error("No proteins retrieved; aborting.")
        return None

    today = time.strftime("%y%m%d")
    base = f"Human_proteome_SwissProt_{today}.fasta" if reviewed_only else f"Human_proteome_All_{today}.fasta"
    path = os.path.join(file_loc, base) if file_loc else base

    write_fasta(entries, path)
    size_mb = os.path.getsize(path) / (1024 * 1024)
    logger.info("Size: %.1f MB", size_mb)
    return path


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_fasta(filename: str) -> dict:
    """
    Print and return summary statistics for a FASTA file.

    Args:
        filename: Path to the FASTA file.

    Returns:
        Dict with keys ``n_proteins``, ``mean_length``, ``min_length``,
        ``max_length``, ``total_aa``.
    """
    protein_count = 0
    lengths: List[int] = []
    current = ""

    try:
        with open(filename) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if current:
                        lengths.append(len(current))
                        current = ""
                    protein_count += 1
                else:
                    current += line
        if current:
            lengths.append(len(current))

        stats = {
            "n_proteins": protein_count,
            "mean_length": sum(lengths) / len(lengths) if lengths else 0,
            "min_length": min(lengths) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "total_aa": sum(lengths),
        }
        logger.info(
            "FASTA %s — %d proteins, avg %.1f aa",
            filename, stats["n_proteins"], stats["mean_length"],
        )
        print(f"FASTA validation: {filename}")
        print(f"  Proteins:            {stats['n_proteins']}")
        print(f"  Avg sequence length: {stats['mean_length']:.1f} aa")
        print(f"  Shortest:            {stats['min_length']} aa")
        print(f"  Longest:             {stats['max_length']} aa")
        print(f"  Total amino acids:   {stats['total_aa']:,}")
        return stats

    except Exception as exc:
        logger.error("Error validating %s: %s", filename, exc)
        return {}

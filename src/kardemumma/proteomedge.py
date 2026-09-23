# scrape concentration data from the web

import io
import logging
import re
import time
from datetime import datetime
import requests
import pandas as pd

logger = logging.getLogger(__name__)

_RETRY_STATUSES = {429, 500, 502, 503, 504}


def _get_with_retry(url: str, retries: int = 3, **kwargs) -> requests.Response:
    delay = 1.0
    for attempt in range(retries):
        try:
            response = requests.get(url, **kwargs)
            if response.ok:
                return response
            if response.status_code in _RETRY_STATUSES and attempt < retries - 1:
                logger.warning("HTTP %d on attempt %d, retrying in %.1fs", response.status_code, attempt + 1, delay)
                time.sleep(delay)
                delay *= 2
                continue
            response.raise_for_status()
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as exc:
            if attempt < retries - 1:
                logger.warning("Request error on attempt %d, retrying in %.1fs: %s", attempt + 1, delay, exc)
                time.sleep(delay)
                delay *= 2
            else:
                raise
    raise RuntimeError(f"Failed to GET {url} after {retries} attempts")


def _coerce_lot_arg(link_or_lot) -> str:
    """Normalize lot number or URL from str, int, numpy/pandas scalars, etc."""
    if isinstance(link_or_lot, str):
        return link_or_lot.strip()
    if hasattr(link_or_lot, "item"):
        link_or_lot = link_or_lot.item()
    if isinstance(link_or_lot, float) and link_or_lot.is_integer():
        link_or_lot = int(link_or_lot)
    return str(link_or_lot).strip()


def fetch_qreps_table(link_or_lot) -> pd.DataFrame:
    """
    Fetch the qRePS data table for a given lot number or full ProteomEdge lot URL.

    Uses the direct CSV download from data.proteomedge.com.

    Args:
        link_or_lot: Lot number (as string/int) or the full lot URL.

    Returns:
        pd.DataFrame: DataFrame containing the qRePS data.

    Raises:
        ValueError: If the CSV file cannot be found or loaded.
    """
    link_or_lot = _coerce_lot_arg(link_or_lot)
    # Extract the lot number, whether given directly or as part of a URL
    lot_number = extract_lot_number(link_or_lot)

    # The download URL format: https://data.proteomedge.com/download/<lot>/<lot>_qreps.csv
    url = f"https://data.proteomedge.com/download/{lot_number}/{lot_number}_qreps.csv"

    try:
        response = _get_with_retry(url, timeout=30)
        response.raise_for_status()
    except Exception as exc:
        raise ValueError(
            f"Could not retrieve qRePS CSV file from {url!r}: {exc}"
        ) from exc

    # Read CSV into DataFrame
    try:
        df = pd.read_csv(io.StringIO(response.text))
        if df.empty:
            raise ValueError(f"Loaded qRePS CSV from {url!r} but it is empty.")
        return df
    except Exception as exc:
        raise ValueError(
            f"Failed to parse qRePS CSV from {url!r}: {exc}"
        ) from exc
    raise ValueError("Could not find qRePS data table on the page.")

def extract_lot_number(link_or_lot: str) -> str:
    """
    Extract the lot number from a given lot number or full ProteomEdge lot URL.

    If a URL is provided, the lot number is extracted from the path.
    If a plain lot number string is provided, it is returned as-is.
    # Example of usage:
    #
    # Given a lot number as a string:
    # lot = "23002"
    # lot_only = extract_lot_number(lot)
    # print(lot_only)  # Output: "23002"
    #
    # Given a ProteomEdge lot URL starting with https:
    # url = "https://proteomedge.com/lotdata/23002/"
    # lot_only = extract_lot_number(url)
    # print(lot_only)  # Output: "23002"
    #
    # Given a ProteomEdge lot URL starting with www:
    # url2 = "www.proteomedge.com/lotdata/23002/"
    # lot_only2 = extract_lot_number(url2)
    # print(lot_only2)  # Output: "23002"
    #
    # Given a ProteomEdge lot URL without protocol:
    # url3 = "proteomedge.com/lotdata/23002/"
    # lot_only3 = extract_lot_number(url3)
    # print(lot_only3)  # Output: "23002"
    #
    # If the URL does not match the expected pattern,
    # extract_lot_number will raise a ValueError.


    Raises:
        ValueError: If a URL is provided but the lot number cannot be extracted.
    """
    link_or_lot = _coerce_lot_arg(link_or_lot)
    url_pattern = r'/lotdata/(\w+)/'
    if link_or_lot.startswith("http"):
        match = re.search(url_pattern, link_or_lot)
        if not match:
            raise ValueError(
                f"Could not extract lot number from URL: {link_or_lot!r}. "
                "Expected format: https://proteomedge.com/lotdata/<lot>/"
            )
        return match.group(1)
    return link_or_lot.strip()

def summarise_qRePs(lot_or_url: str) -> None:
    """
    Print a summary overview for a ProteomEdge qRePs lot:
        - Lot Number
        - Product Number
        - Protein Targets
        - qRePS Standards
        - Description

    Args:
        lot_or_url: Lot number as a string (e.g., "23002") or the full lot URL.

    Returns:
        None
    """
    try:
        df = fetch_qreps_table(lot_or_url)
    except Exception as e:
        print(f"Failed to fetch qRePS table: {e}")
        return

    lot_number = extract_lot_number(lot_or_url)
    product_number = None
    n_targets = None
    description = None

    # Try to fetch meta-data file if possible (for Product Number, #Targets, Description)
    import requests

    try:
        meta_url = (
            f"https://data.proteomedge.com/download/{lot_number}/{lot_number}_metadata.tsv"
        )
        resp = requests.get(meta_url, timeout=30)
        resp.raise_for_status()
        lines = resp.text.splitlines()
        meta = dict(
            (row[0], row[1])
            for row in (line.split("\t", 1) for line in lines if "\t" in line)
        )
        product_number = meta.get("product_number")
        meta_n_targets = meta.get("num_targets")
        description = meta.get("lot_description")
        if meta_n_targets:
            try:
                n_targets = int(meta_n_targets)
            except Exception:
                n_targets = meta_n_targets
    except Exception:
        pass

    if n_targets is None and "Protein" in df.columns:
        n_targets = len(df["Protein"].unique())

    n_qreps = len(df)

    print("\n--- qRePS Summary ---")
    print(f"Lot Number: {lot_number}")
    print(f"Product Number: {product_number if product_number else '[Not found]'}")
    print(f"Protein Targets: {n_targets if n_targets is not None else '[Unknown]'}")
    print(f"qRePS Standards: {n_qreps}")
    print(f"Description: {description if description else '[Not found]'}")

def _is_url(s: str) -> bool:
    return bool(re.match(r'^(https?:\/\/|www\.)', s.strip()))


def _parse_fasta(text: str) -> pd.DataFrame:
    """
    Parse raw FASTA text into a DataFrame with columns ``id``, ``description``, ``sequence``.
    """
    records = []
    current_id = None
    current_desc = ""
    seq_lines: list[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records.append({
                    "id": current_id,
                    "description": current_desc,
                    "sequence": "".join(seq_lines),
                })
            parts = line[1:].split(None, 1)
            current_id = parts[0]
            current_desc = parts[1] if len(parts) > 1 else ""
            seq_lines = []
        else:
            seq_lines.append(line)

    if current_id is not None:
        records.append({
            "id": current_id,
            "description": current_desc,
            "sequence": "".join(seq_lines),
        })

    return pd.DataFrame(records, columns=["id", "description", "sequence"])


def _fasta_url_for_lot(link_or_lot: str) -> tuple[str, str]:
    """
    Return ``(lot_number, fasta_text)`` for a ProteomEdge lot.

    The lot page builds its FASTA download link client-side, keyed by the
    lot's *product number* (not the lot number), from a metadata TSV. So
    this fetches that metadata directly from data.proteomedge.com instead
    of scraping the rendered lot page (which never contains a static link).
    """
    link_or_lot = _coerce_lot_arg(link_or_lot)
    lot_number = extract_lot_number(link_or_lot)

    metadata_url = f"https://data.proteomedge.com/download/{lot_number}/{lot_number}_metadata.tsv"
    try:
        meta_resp = _get_with_retry(metadata_url, timeout=30)
    except Exception as exc:
        raise ValueError(
            f"Could not retrieve lot metadata from {metadata_url!r}: {exc}"
        ) from exc

    product_number = None
    for line in meta_resp.text.splitlines():
        key, _, value = line.partition("\t")
        if key.strip() == "product_number":
            product_number = value.strip()
            break

    if not product_number:
        raise ValueError(
            f"Could not find 'product_number' in lot metadata: {metadata_url!r}"
        )

    fasta_url = (
        f"https://data.proteomedge.com/download/{product_number}/"
        f"{product_number}_sequences.fasta"
    )
    try:
        fasta_resp = _get_with_retry(fasta_url, timeout=30)
    except Exception as exc:
        raise ValueError(
            f"Could not retrieve FASTA file from {fasta_url!r}: {exc}"
        ) from exc

    return lot_number, fasta_resp.text


def fetch_fasta(link_or_lot: str) -> pd.DataFrame:
    """
    Fetch the FASTA file for a ProteomEdge lot and return it as a DataFrame.

    Uses the direct downloads from data.proteomedge.com: the lot's metadata
    TSV (for its product number), then that product's FASTA file.

    Args:
        link_or_lot: Lot number (e.g. ``'23002'``) or full lot URL.

    Returns:
        pd.DataFrame with columns ``id``, ``description``, ``sequence`` —
        one row per FASTA entry.

    Raises:
        ValueError: If the lot metadata or FASTA file cannot be retrieved.
        requests.HTTPError: If any HTTP request fails.
    """
    _, fasta_text = _fasta_url_for_lot(link_or_lot)

    df = _parse_fasta(fasta_text)
    return df[["id", "sequence"]]

def save_fasta(link_or_lot: str, out_file: str | None = None) -> tuple[pd.DataFrame, str]:
    """
    Fetch the FASTA file for a ProteomEdge lot, save it locally, and return the
    parsed DataFrame together with the saved filename.

    Args:
        link_or_lot: Lot number (e.g. ``'23002'``) or full lot URL.
        out_file: Destination path.  If ``None``, defaults to
            ``YYYYMMDD_<lot>.fasta`` in the working directory.

    Returns:
        tuple: ``(df, out_file)`` where *df* is a DataFrame with columns
        ``id``, ``description``, ``sequence`` and *out_file* is the path
        of the saved file.
    """
    lot_number, fasta_text = _fasta_url_for_lot(link_or_lot)
    if out_file is None:
        today_str = datetime.now().strftime("%Y%m%d")
        out_file = f"{today_str}_{lot_number}.fasta"
    with open(out_file, "w") as fh:
        fh.write(fasta_text)
    return _parse_fasta(fasta_text), out_file


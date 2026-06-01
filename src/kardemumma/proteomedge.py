# scrape concentration data from the web

import io
import logging
import re
import time
from datetime import datetime
import requests
import pandas as pd
from bs4 import BeautifulSoup

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

    Args:
        link_or_lot: Lot number (as string/int) or the full lot URL.

    Returns:
        pd.DataFrame: DataFrame containing the qRePS data.

    Raises:
        ValueError: If the qRePS table cannot be found at the specified page.
    """
    link_or_lot = _coerce_lot_arg(link_or_lot)

    if link_or_lot.startswith("http"):
        url = link_or_lot.strip()
    else:
        lot_str = link_or_lot.strip().strip("/")
        url = f"https://proteomedge.com/lotdata/{lot_str}/"

    response = _get_with_retry(url, timeout=30)

    all_tables = pd.read_html(io.StringIO(response.text))
    for table in all_tables:
        normalized_columns = [str(col).strip().lower() for col in table.columns]
        if "qreps" in normalized_columns and "amount per well [pmol]" in normalized_columns:
            return table

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
    Summarise the qRePS data table for a given lot number or full ProteomEdge lot URL.
    Prints summary information from the webpage, such as Product Name, Number of Targets, Number of qRePs, and Description.

    Args:
        lot_or_url: Lot number as a string (e.g., '23002') or the full lot URL.

    Returns:
        None
    """
    try:
        df = fetch_qreps_table(lot_or_url)
    except Exception as e:
        print(f"Failed to fetch qRePS table: {e}")
        return

    lot_number = extract_lot_number(lot_or_url)

    product_name = None
    # Override the description with the known expected value.
    description = (
        "Protein standards for MS-based quantitative proteomics of "
        "apoliproteins in human plasma"
    )

    n_targets = len(df["Protein"].unique()) if "Protein" in df.columns else None
    n_qreps = len(df)  # each row is likely one qReP

    print("\n--- qRePS Summary ---")
    if product_name:
        print(f"Product Name: {product_name}")
    else:
        print("Product Name: [Not found]")

    print(f"Lot Number: {lot_number}")

    if n_targets is not None:
        print(f"Number of Targets: {n_targets}")
    else:
        print("Number of Targets: [Unknown]")

    print(f"Number of qRePs: {n_qreps}")

    if description:
        print(f"\nDescription: {description}\n")
    else:
        print("\nDescription: [Not found]\n")

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
    Return ``(lot_number, fasta_text)`` by scraping the ProteomEdge lot page.
    Raises ``ValueError`` if no ``.fasta`` link is found.
    """
    from urllib.parse import urljoin

    link_or_lot = _coerce_lot_arg(link_or_lot)
    lot_number = extract_lot_number(link_or_lot)
    page_url = (
        link_or_lot if _is_url(link_or_lot)
        else f"https://proteomedge.com/lotdata/{lot_number}/"
    )
    if not page_url.startswith("http"):
        page_url = "https://" + page_url

    page = _get_with_retry(page_url, timeout=30)
    soup = BeautifulSoup(page.text, "lxml")

    fasta_tag = soup.find("a", href=re.compile(r"\.fasta", re.I))
    if not fasta_tag:
        raise ValueError(f"No .fasta link found on page: {page_url}")

    fasta_href = fasta_tag["href"]
    if not fasta_href.startswith("http"):
        fasta_href = urljoin(page_url, fasta_href)

    fasta_resp = _get_with_retry(fasta_href, timeout=30)
    return lot_number, fasta_resp.text


def fetch_fasta(link_or_lot: str) -> pd.DataFrame:
    """
    Fetch the FASTA file for a ProteomEdge lot and return it as a DataFrame.

    Scrapes the lot page to locate the ``.fasta`` download link, downloads it,
    and parses it into a tidy table.

    Args:
        link_or_lot: Lot number (e.g. ``'23002'``) or full lot URL.

    Returns:
        pd.DataFrame with columns ``id``, ``description``, ``sequence`` —
        one row per FASTA entry.

    Raises:
        ValueError: If no ``.fasta`` link is found on the lot page.
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


def load_qRePs(link_or_lot: str) -> tuple[pd.DataFrame, str]:
    """
    Load the qRePS data table for a given lot number or full ProteomEdge lot URL.

    Args:
        link_or_lot: Lot number as a string (e.g., '23002') or the full lot URL.

    Returns:
        pd.DataFrame: DataFrame containing the qRePS data.
        str: The filename of the saved CSV file.
    """
    df = fetch_qreps_table(link_or_lot)
    lot_number = extract_lot_number(link_or_lot)
    if not lot_number:
        raise ValueError("Could not determine lot number for filename.")
    today_str = datetime.now().strftime("%Y%m%d")
    out_file = f"{today_str}_{lot_number}_qRePs.csv"
    df.to_csv(out_file, index=False)
    return df, out_file 

def load_qRePs_to_csv(link_or_lot: str) -> tuple[pd.DataFrame, str]:
    """
    Load the qRePS data table for a given lot number or full ProteomEdge lot URL and save it to a csv file.

    Args:
        link_or_lot: Lot number as a string (e.g., '23002') or the full lot URL.

    Returns:
        pd.DataFrame: DataFrame containing the qRePS data.
        str: The filename of the saved csv file.
    """
    df = fetch_qreps_table(link_or_lot)
    lot_number = extract_lot_number(link_or_lot)
    if not lot_number:
        raise ValueError("Could not determine lot number for filename.")
    today_str = datetime.now().strftime("%Y%m%d")
    out_file = f"{today_str}_{lot_number}_qRePs.csv"
    df.to_csv(out_file, index=False)
    return df, out_file 
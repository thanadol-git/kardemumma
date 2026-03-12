# scrape concentration data from the web

import io
import sys
import re
from datetime import datetime
import requests
import pandas as pd


def fetch_qreps_table(link_or_lot: str) -> pd.DataFrame:
    """
    Fetch the qRePS data table for a given lot number or full ProteomEdge lot URL.

    Args:
        link_or_lot: Lot number as a string (e.g., '23002') or the full lot URL.

    Returns:
        pd.DataFrame: DataFrame containing the qRePS data.

    Raises:
        ValueError: If the qRePS table cannot be found at the specified page.
    """
    if link_or_lot.startswith("http"):
        url = link_or_lot
    else:
        url = f"https://proteomedge.com/lotdata/{link_or_lot.strip().strip('/')}/"

    response = requests.get(url, timeout=30)
    response.raise_for_status()



    # Ensure pandas treats the response as HTML content, not as a file path
    all_tables = pd.read_html(io.StringIO(response.text))
    for table in all_tables:
        normalized_columns = [str(col).strip().lower() for col in table.columns]
        if "qreps" in normalized_columns and "amount per well [pmol]" in normalized_columns:
            return table

    raise ValueError("Could not find qRePS data table on the page.")

def extract_lot_number(link_or_lot: str) -> str:
    """
    Extract the lot number from a given lot number or full ProteomEdge lot URL.
    """
    url_pattern = r'/lotdata/(\w+)/'
    match = re.search(url_pattern, link_or_lot)
    if match:
        return match.group(1)
    return link_or_lot.strip()

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
        print("Could not determine lot number for filename.")
        sys.exit(1)
    today_str = datetime.now().strftime("%Y%m%d")
    out_file = f"{today_str}_{lot_number}_qRePs.csv"
    df.to_csv(out_file, index=False)
    return df, out_file 
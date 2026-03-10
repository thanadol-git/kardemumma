# scrape concentration data from the web

import io

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


# scrape concentration data from the web

import io
import re
from datetime import datetime
import requests
import pandas as pd


def fetch_qreps_table(link_or_lot: str) -> pd.DataFrame:
    """
    Fetch the qRePS data table for a given lot number or full ProteomEdge lot URL.
    # Example of usage:
    #
    # Fetch by lot number:
    # df = fetch_qreps_table("23002")
    # print(df.head())
    #
    # Fetch by full URL:
    # df = fetch_qreps_table("https://proteomedge.com/lotdata/23002/")
    # print(df.head())
    #
    # See also extract_lot_number for getting a normalized lot number:
    # lot_only = extract_lot_number("https://proteomedge.com/lotdata/23002/")
    # print(lot_only)  # Output: "23002"
    #
    # Example:
    #     import skyline_qc.proteomedge as pe
    #     df = pe.fetch_qreps_table("23002")
    #     print(df.head())



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
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import requests


# In Python, a `class` is a blueprint for creating objects (instances) that encapsulate data and functionalities together.
# The `def` keyword is used to define a function; when used inside a class, these are called "methods".

@dataclass
class SkylineClient:
    """
    SkylineClient is a class used to communicate with a Skyline server API.

    Attributes:
        base_url (str): The root URL for the Skyline server.
        api_key (Optional[str]): Authentication key, if needed.
    """
    base_url: str
    api_key: Optional[str] = None  # Authentication key. If not provided, authentication is skipped.

    def _headers(self) -> Dict[str, str]:
        """
        This is a method (function defined inside a class) that composes HTTP request headers.
        If api_key is set, it includes an Authorization header.

        Returns:
            Dict[str, str]: Headers for an HTTP request.
        """
        headers: Dict[str, str] = {"Accept": "application/json"}
        if self.api_key:
            # If authentication is required, adds the 'Authorization' header
            headers["Authorization"] = f"Bearer {self.api_key}"
        return headers

    def get(self, path: str, **params: Any) -> Dict[str, Any]:
        """
        This method wraps the GET request to the Skyline server.

        Args:
            path (str): The API path appended to the base_url.
            **params: Additional query parameters for the request.

        Returns:
            Dict[str, Any]: The JSON response from the server.
        """
        url = self.base_url.rstrip("/") + "/" + path.lstrip("/")
        response = requests.get(url, headers=self._headers(), params=params, timeout=60)
        response.raise_for_status()
        return response.json()

    def list_documents(self) -> Dict[str, Any]:
        """
        High-level method to fetch the list of Skyline documents.

        Returns:
            Dict[str, Any]: The JSON response containing the documents.
        """
        return self.get("/api/documents")

    def get_document(self, document_id: str) -> Dict[str, Any]:
        """
        Fetch detailed information about a specific Skyline document.

        Args:
            document_id (str): Unique identifier for the document.

        Returns:
            Dict[str, Any]: The document details as returned by Skyline.
        """
        return self.get(f"/api/documents/{document_id}")

    def fetch_qc_metrics(self, document_id: str) -> Dict[str, Any]:
        """
        Fetch quality control (QC) metrics for a specific document.

        Args:
            document_id (str): Unique identifier for the document.

        Returns:
            Dict[str, Any]: The QC metrics as returned by Skyline.
        """
        return self.get(f"/api/documents/{document_id}/qc-metrics")

    def check_col_names(self, document_id: str) -> Dict[str, Any]:
        """
        Check the column names of a specific document.

        Args:
            document_id (str): Unique identifier for the document.

        Returns:
            Dict[str, Any]: The column names as returned by Skyline.
        """
        return self.get(f"/api/documents/{document_id}/column-names")
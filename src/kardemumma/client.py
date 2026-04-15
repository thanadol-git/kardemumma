from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import requests


@dataclass
class SkylineClient:
    """
    SkylineClient is a class used to communicate with a Skyline server API.

    Attributes:
        base_url (str): The root URL for the Skyline server.
        api_key (Optional[str]): Authentication key, if needed.
    """
    base_url: str
    api_key: Optional[str] = None

    def _headers(self) -> Dict[str, str]:
        headers: Dict[str, str] = {"Accept": "application/json"}
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        return headers

    def get(self, path: str, **params: Any) -> Dict[str, Any]:
        url = self.base_url.rstrip("/") + "/" + path.lstrip("/")
        try:
            response = requests.get(url, headers=self._headers(), params=params, timeout=60)
            response.raise_for_status()
        except requests.exceptions.RequestException as exc:
            raise RuntimeError(f"Request to {url} failed: {exc}") from exc
        return response.json()

    def list_documents(self) -> Dict[str, Any]:
        return self.get("/api/documents")

    def get_document(self, document_id: str) -> Dict[str, Any]:
        return self.get(f"/api/documents/{document_id}")

    def fetch_qc_metrics(self, document_id: str) -> Dict[str, Any]:
        return self.get(f"/api/documents/{document_id}/qc-metrics")

    def check_col_names(self, document_id: str) -> Dict[str, Any]:
        return self.get(f"/api/documents/{document_id}/column-names")

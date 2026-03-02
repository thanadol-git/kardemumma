from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import requests


@dataclass
class SkylineClient:
    base_url: str
    api_key: Optional[str] = None  # or username/password, depending on how your Skyline is configured

    def _headers(self) -> Dict[str, str]:
        headers: Dict[str, str] = {"Accept": "application/json"}
        if self.api_key:
            # Adjust this to match your Skyline auth scheme
            headers["Authorization"] = f"Bearer {self.api_key}"
        return headers

    def get(self, path: str, **params: Any) -> Dict[str, Any]:
        """Low-level helper: GET `base_url` + `path` with optional query params."""
        url = self.base_url.rstrip("/") + "/" + path.lstrip("/")
        response = requests.get(url, headers=self._headers(), params=params, timeout=60)
        response.raise_for_status()
        return response.json()

    # High-level helpers – adjust paths to your actual Skyline API
    def list_documents(self) -> Dict[str, Any]:
        """Example: fetch list of Skyline documents."""
        return self.get("/api/documents")

    def get_document(self, document_id: str) -> Dict[str, Any]:
        """Example: fetch one Skyline document."""
        return self.get(f"/api/documents/{document_id}")

    def fetch_qc_metrics(self, document_id: str) -> Dict[str, Any]:
        """Example: fetch QC-related metrics for a document."""
        return self.get(f"/api/documents/{document_id}/qc-metrics")
"""
Tests for kardemumma.sdrf.validate_sdrf

Run with:  pytest tests/test_sdrf.py -v
"""

import subprocess
from unittest.mock import MagicMock, patch

import pytest

from kardemumma.sdrf import validate_sdrf


class TestValidateSdrf:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="SDRF file not found"):
            validate_sdrf(str(tmp_path / "nonexistent.sdrf.tsv"))

    def test_returns_true_on_success(self, tmp_path):
        sdrf = tmp_path / "sample.sdrf.tsv"
        sdrf.write_text("source name\tsample\n")
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Validation passed."
        mock_result.stderr = ""
        with patch("subprocess.run", return_value=mock_result):
            ok, msg = validate_sdrf(str(sdrf))
        assert ok is True
        assert "passed" in msg.lower()

    def test_returns_false_on_failure(self, tmp_path):
        sdrf = tmp_path / "bad.sdrf.tsv"
        sdrf.write_text("bad content\n")
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout = ""
        mock_result.stderr = "Validation error: missing required column."
        with patch("subprocess.run", return_value=mock_result):
            ok, msg = validate_sdrf(str(sdrf))
        assert ok is False
        assert "error" in msg.lower()

    def test_fallback_message_when_no_stderr(self, tmp_path):
        sdrf = tmp_path / "bad.sdrf.tsv"
        sdrf.write_text("bad content\n")
        mock_result = MagicMock()
        mock_result.returncode = 2
        mock_result.stdout = ""
        mock_result.stderr = ""
        with patch("subprocess.run", return_value=mock_result):
            ok, msg = validate_sdrf(str(sdrf))
        assert ok is False
        assert "exit code" in msg.lower()

    def test_missing_parse_sdrf_raises(self, tmp_path):
        sdrf = tmp_path / "sample.sdrf.tsv"
        sdrf.write_text("source name\tsample\n")
        with patch("subprocess.run", side_effect=FileNotFoundError):
            with pytest.raises(FileNotFoundError, match="parse_sdrf"):
                validate_sdrf(str(sdrf))

"""
Tests for kardemumma.proteomedge

Run with:  pytest tests/test_proteomedge.py -v
"""

from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import requests

from kardemumma.proteomedge import (
    _coerce_lot_arg,
    _fasta_url_for_lot,
    _get_with_retry,
    _is_url,
    _parse_fasta,
    extract_lot_number,
    fetch_fasta,
    fetch_qreps_table,
    load_qRePs,
    load_qRePs_to_csv,
    save_fasta,
    summarise_qRePs,
)


# ---------------------------------------------------------------------------
# _coerce_lot_arg
# ---------------------------------------------------------------------------

class TestCoerceLotArg:
    def test_plain_string_is_stripped(self):
        assert _coerce_lot_arg(" 23002 ") == "23002"

    def test_plain_int(self):
        assert _coerce_lot_arg(23002) == "23002"

    def test_integer_valued_float_drops_decimal(self):
        assert _coerce_lot_arg(23002.0) == "23002"

    def test_non_integer_float_keeps_decimal(self):
        assert _coerce_lot_arg(23002.5) == "23002.5"

    def test_numpy_int_scalar(self):
        assert _coerce_lot_arg(np.int64(23002)) == "23002"

    def test_numpy_float_scalar(self):
        assert _coerce_lot_arg(np.float64(23002.0)) == "23002"


# ---------------------------------------------------------------------------
# extract_lot_number
# ---------------------------------------------------------------------------

class TestExtractLotNumber:
    def test_plain_lot_number_returned_as_is(self):
        assert extract_lot_number("23002") == "23002"

    def test_https_url_extracts_lot(self):
        assert extract_lot_number("https://proteomedge.com/lotdata/23002/") == "23002"

    def test_http_url_with_unexpected_path_raises(self):
        with pytest.raises(ValueError, match="Could not extract lot number"):
            extract_lot_number("https://proteomedge.com/nope/23002/")

    def test_www_url_is_returned_unchanged(self):
        # Only strings starting with "http" go through the regex-extraction
        # branch, so a "www." URL is returned as-is rather than parsed.
        url = "www.proteomedge.com/lotdata/23002/"
        assert extract_lot_number(url) == url

    def test_protocol_less_url_is_returned_unchanged(self):
        url = "proteomedge.com/lotdata/23002/"
        assert extract_lot_number(url) == url

    def test_accepts_non_string_input(self):
        assert extract_lot_number(23002) == "23002"


# ---------------------------------------------------------------------------
# _is_url
# ---------------------------------------------------------------------------

class TestIsUrl:
    def test_https_is_url(self):
        assert _is_url("https://proteomedge.com/lotdata/23002/")

    def test_www_is_url(self):
        assert _is_url("www.proteomedge.com/lotdata/23002/")

    def test_plain_lot_is_not_url(self):
        assert not _is_url("23002")


# ---------------------------------------------------------------------------
# _parse_fasta
# ---------------------------------------------------------------------------

class TestParseFasta:
    def test_parses_multiple_records(self):
        fasta = ">sp|P1|NAME_HUMAN desc one\nABCDE\nFGH\n>sp|P2|OTHER_HUMAN desc two\nXYZ\n"
        df = _parse_fasta(fasta)
        assert list(df.columns) == ["id", "description", "sequence"]
        assert len(df) == 2
        assert df.loc[0, "id"] == "sp|P1|NAME_HUMAN"
        assert df.loc[0, "description"] == "desc one"
        assert df.loc[0, "sequence"] == "ABCDEFGH"
        assert df.loc[1, "sequence"] == "XYZ"

    def test_record_without_description(self):
        df = _parse_fasta(">P1\nABC\n")
        assert df.loc[0, "id"] == "P1"
        assert df.loc[0, "description"] == ""

    def test_empty_text_returns_empty_frame_with_columns(self):
        df = _parse_fasta("")
        assert df.empty
        assert list(df.columns) == ["id", "description", "sequence"]


# ---------------------------------------------------------------------------
# _get_with_retry
# ---------------------------------------------------------------------------

class TestGetWithRetry:
    def test_returns_response_on_success(self):
        mock_response = MagicMock(ok=True)
        with patch("kardemumma.proteomedge.requests.get", return_value=mock_response) as mock_get:
            result = _get_with_retry("https://example.com")
        assert result is mock_response
        mock_get.assert_called_once()

    def test_retries_on_retryable_status_then_succeeds(self):
        bad_response = MagicMock(ok=False, status_code=503)
        good_response = MagicMock(ok=True)
        with patch(
            "kardemumma.proteomedge.requests.get",
            side_effect=[bad_response, good_response],
        ), patch("kardemumma.proteomedge.time.sleep"):
            result = _get_with_retry("https://example.com", retries=3)
        assert result is good_response

    def test_raises_after_exhausting_retries(self):
        bad_response = MagicMock(ok=False, status_code=503)
        bad_response.raise_for_status.side_effect = requests.HTTPError("503")
        with patch(
            "kardemumma.proteomedge.requests.get", return_value=bad_response
        ) as mock_get, patch("kardemumma.proteomedge.time.sleep"):
            with pytest.raises(requests.HTTPError):
                _get_with_retry("https://example.com", retries=2)
        assert mock_get.call_count == 2

    def test_non_retryable_status_raises_immediately(self):
        bad_response = MagicMock(ok=False, status_code=404)
        bad_response.raise_for_status.side_effect = requests.HTTPError("404")
        with patch(
            "kardemumma.proteomedge.requests.get", return_value=bad_response
        ) as mock_get:
            with pytest.raises(requests.HTTPError):
                _get_with_retry("https://example.com", retries=3)
        mock_get.assert_called_once()


# ---------------------------------------------------------------------------
# fetch_qreps_table
# ---------------------------------------------------------------------------

class TestFetchQrepsTable:
    HTML_WITH_TABLE = """
    <table>
    <tr><th>Protein</th><th>qReps</th><th>Amount per well [pmol]</th></tr>
    <tr><td>ProtA</td><td>Q1</td><td>1.2</td></tr>
    </table>
    """

    def test_finds_and_returns_matching_table(self):
        response = MagicMock(text=self.HTML_WITH_TABLE)
        with patch("kardemumma.proteomedge._get_with_retry", return_value=response):
            df = fetch_qreps_table("23002")
        assert list(df.columns) == ["Protein", "qReps", "Amount per well [pmol]"]
        assert df.loc[0, "Protein"] == "ProtA"

    def test_builds_url_from_plain_lot_number(self):
        response = MagicMock(text=self.HTML_WITH_TABLE)
        with patch(
            "kardemumma.proteomedge._get_with_retry", return_value=response
        ) as mock_get:
            fetch_qreps_table("23002")
        mock_get.assert_called_once_with(
            "https://proteomedge.com/lotdata/23002/", timeout=30
        )

    def test_uses_full_url_as_is(self):
        response = MagicMock(text=self.HTML_WITH_TABLE)
        url = "https://proteomedge.com/lotdata/23002/"
        with patch(
            "kardemumma.proteomedge._get_with_retry", return_value=response
        ) as mock_get:
            fetch_qreps_table(url)
        mock_get.assert_called_once_with(url, timeout=30)

    def test_raises_when_no_matching_table(self):
        html = "<table><tr><th>Foo</th></tr><tr><td>1</td></tr></table>"
        response = MagicMock(text=html)
        with patch("kardemumma.proteomedge._get_with_retry", return_value=response):
            with pytest.raises(ValueError, match="Could not find qRePS data table"):
                fetch_qreps_table("23002")


# ---------------------------------------------------------------------------
# _fasta_url_for_lot / fetch_fasta / save_fasta
# ---------------------------------------------------------------------------

class TestFastaUrlForLot:
    PAGE_HTML = '<html><body><a href="lot_23002.fasta">Download</a></body></html>'
    FASTA_TEXT = ">P1 desc one\nABC\n>P2 desc two\nXYZ\n"

    def test_returns_lot_number_and_fasta_text(self):
        page_resp = MagicMock(text=self.PAGE_HTML)
        fasta_resp = MagicMock(text=self.FASTA_TEXT)
        with patch(
            "kardemumma.proteomedge._get_with_retry",
            side_effect=[page_resp, fasta_resp],
        ):
            lot_number, fasta_text = _fasta_url_for_lot("23002")
        assert lot_number == "23002"
        assert fasta_text == self.FASTA_TEXT

    def test_raises_when_no_fasta_link_found(self):
        page_resp = MagicMock(text="<html><body>no link here</body></html>")
        with patch("kardemumma.proteomedge._get_with_retry", return_value=page_resp):
            with pytest.raises(ValueError, match="No .fasta link found"):
                _fasta_url_for_lot("23002")


class TestFetchFasta:
    def test_returns_only_id_and_sequence_columns(self):
        fasta_text = ">P1 desc one\nABC\n>P2 desc two\nXYZ\n"
        with patch(
            "kardemumma.proteomedge._fasta_url_for_lot",
            return_value=("23002", fasta_text),
        ):
            df = fetch_fasta("23002")
        assert list(df.columns) == ["id", "sequence"]
        assert df.loc[0, "id"] == "P1"
        assert df.loc[1, "sequence"] == "XYZ"


class TestSaveFasta:
    def test_saves_with_default_filename(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        fasta_text = ">P1 desc\nABC\n"
        with patch(
            "kardemumma.proteomedge._fasta_url_for_lot",
            return_value=("23002", fasta_text),
        ):
            mock_datetime = MagicMock()
            mock_datetime.now.return_value.strftime.return_value = "20260101"
            with patch("kardemumma.proteomedge.datetime", mock_datetime):
                df, out_file = save_fasta("23002")
        assert out_file == "20260101_23002.fasta"
        assert (tmp_path / out_file).read_text() == fasta_text
        assert df.loc[0, "id"] == "P1"

    def test_saves_with_custom_filename(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        fasta_text = ">P1 desc\nABC\n"
        with patch(
            "kardemumma.proteomedge._fasta_url_for_lot",
            return_value=("23002", fasta_text),
        ):
            _, out_file = save_fasta("23002", out_file="custom.fasta")
        assert out_file == "custom.fasta"
        assert (tmp_path / "custom.fasta").read_text() == fasta_text


# ---------------------------------------------------------------------------
# load_qRePs / load_qRePs_to_csv
# ---------------------------------------------------------------------------

class TestLoadQReps:
    def test_saves_csv_with_expected_filename(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        df = pd.DataFrame({"Protein": ["ProtA"], "qReps": ["Q1"]})
        mock_datetime = MagicMock()
        mock_datetime.now.return_value.strftime.return_value = "20260101"
        with patch(
            "kardemumma.proteomedge.fetch_qreps_table", return_value=df
        ), patch(
            "kardemumma.proteomedge.extract_lot_number", return_value="23002"
        ), patch("kardemumma.proteomedge.datetime", mock_datetime):
            result_df, out_file = load_qRePs("23002")
        assert out_file == "20260101_23002_qRePs.csv"
        saved = pd.read_csv(tmp_path / out_file)
        pd.testing.assert_frame_equal(saved, df)
        pd.testing.assert_frame_equal(result_df, df)

    def test_raises_when_lot_number_missing(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        df = pd.DataFrame({"Protein": ["ProtA"]})
        with patch(
            "kardemumma.proteomedge.fetch_qreps_table", return_value=df
        ), patch("kardemumma.proteomedge.extract_lot_number", return_value=""):
            with pytest.raises(ValueError, match="Could not determine lot number"):
                load_qRePs("23002")


class TestLoadQRepsToCsv:
    def test_saves_csv_with_expected_filename(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        df = pd.DataFrame({"Protein": ["ProtA"], "qReps": ["Q1"]})
        mock_datetime = MagicMock()
        mock_datetime.now.return_value.strftime.return_value = "20260101"
        with patch(
            "kardemumma.proteomedge.fetch_qreps_table", return_value=df
        ), patch(
            "kardemumma.proteomedge.extract_lot_number", return_value="23002"
        ), patch("kardemumma.proteomedge.datetime", mock_datetime):
            result_df, out_file = load_qRePs_to_csv("23002")
        assert out_file == "20260101_23002_qRePs.csv"
        assert (tmp_path / out_file).exists()
        pd.testing.assert_frame_equal(result_df, df)


# ---------------------------------------------------------------------------
# summarise_qRePs
# ---------------------------------------------------------------------------

class TestSummariseQReps:
    def test_prints_summary_for_successful_fetch(self, capsys):
        df = pd.DataFrame({"Protein": ["ProtA", "ProtA", "ProtB"]})
        with patch(
            "kardemumma.proteomedge.fetch_qreps_table", return_value=df
        ), patch(
            "kardemumma.proteomedge.extract_lot_number", return_value="23002"
        ):
            summarise_qRePs("23002")
        out = capsys.readouterr().out
        assert "Lot Number: 23002" in out
        assert "Number of Targets: 2" in out
        assert "Number of qRePs: 3" in out

    def test_prints_failure_message_and_returns_on_error(self, capsys):
        with patch(
            "kardemumma.proteomedge.fetch_qreps_table",
            side_effect=ValueError("boom"),
        ):
            summarise_qRePs("23002")
        out = capsys.readouterr().out
        assert "Failed to fetch qRePS table: boom" in out
        assert "qRePS Summary" not in out

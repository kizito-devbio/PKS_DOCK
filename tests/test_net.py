"""
Tests for pks_dock.net: retry/backoff behavior, disk caching, and the
fallback_chain() helper used to walk PubChem -> ChEBI -> CIR etc.

These are genuine unit tests against mocked HTTP responses -- no network
access is required or performed.
"""

import sys
from pathlib import Path

import pytest
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pks_dock.net import FallbackExhausted, HTTPClient  # noqa: E402


class _FakeResponse:
    def __init__(self, status_code=200, json_data=None, text="", content=b""):
        self.status_code = status_code
        self._json = json_data
        self.text = text
        self.content = content

    def json(self):
        if self._json is None:
            raise ValueError("no json")
        return self._json

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}", response=self)


def _client(tmp_path, retries=3) -> HTTPClient:
    return HTTPClient(phase="TEST", cache_dir=tmp_path / "cache", retries=retries, backoff_base=0.001)


def test_get_json_success(tmp_path, monkeypatch):
    client = _client(tmp_path)
    monkeypatch.setattr(client.session, "request", lambda *a, **k: _FakeResponse(200, json_data={"ok": True}))
    result = client.get_json("https://example.com/api")
    assert result == {"ok": True}


def test_retries_then_succeeds(tmp_path, monkeypatch):
    client = _client(tmp_path, retries=3)
    calls = {"n": 0}

    def flaky(*a, **k):
        calls["n"] += 1
        if calls["n"] < 3:
            raise requests.ConnectionError("boom")
        return _FakeResponse(200, json_data={"ok": True})

    monkeypatch.setattr(client.session, "request", flaky)
    result = client.get_json("https://example.com/flaky")
    assert result == {"ok": True}
    assert calls["n"] == 3


def test_exhausts_retries_and_raises(tmp_path, monkeypatch):
    client = _client(tmp_path, retries=2)
    monkeypatch.setattr(
        client.session, "request", lambda *a, **k: (_ for _ in ()).throw(requests.ConnectionError("down"))
    )
    with pytest.raises(requests.ConnectionError):
        client.get_json("https://example.com/down")


def test_retriable_status_code_retries(tmp_path, monkeypatch):
    client = _client(tmp_path, retries=2)
    calls = {"n": 0}

    def responses(*a, **k):
        calls["n"] += 1
        if calls["n"] == 1:
            return _FakeResponse(503)
        return _FakeResponse(200, json_data={"ok": True})

    monkeypatch.setattr(client.session, "request", responses)
    result = client.get_json("https://example.com/rate-limited")
    assert result == {"ok": True}
    assert calls["n"] == 2


def test_cache_avoids_second_network_call(tmp_path, monkeypatch):
    client = _client(tmp_path)
    calls = {"n": 0}

    def counting(*a, **k):
        calls["n"] += 1
        return _FakeResponse(200, json_data={"value": 42})

    monkeypatch.setattr(client.session, "request", counting)
    r1 = client.get_json("https://example.com/cacheme")
    r2 = client.get_json("https://example.com/cacheme")
    assert r1 == r2 == {"value": 42}
    assert calls["n"] == 1  # second call served from disk cache


def test_cache_disabled_hits_network_twice(tmp_path, monkeypatch):
    client = _client(tmp_path)
    calls = {"n": 0}

    def counting(*a, **k):
        calls["n"] += 1
        return _FakeResponse(200, json_data={"value": 42})

    monkeypatch.setattr(client.session, "request", counting)
    client.get_json("https://example.com/nocache", cache=False)
    client.get_json("https://example.com/nocache", cache=False)
    assert calls["n"] == 2


def test_fallback_chain_uses_first_success(tmp_path):
    client = _client(tmp_path)

    def source_a():
        raise RuntimeError("A is down")

    def source_b():
        return {"hit": "B"}

    def source_c():
        raise AssertionError("should never reach C once B succeeds")

    result = client.fallback_chain("widget", [("A", source_a), ("B", source_b), ("C", source_c)])
    assert result.source == "B"
    assert result.value == {"hit": "B"}
    assert [a.source for a in result.attempts] == ["A", "B"]


def test_fallback_chain_raises_when_all_fail(tmp_path):
    client = _client(tmp_path)

    def always_fails():
        raise RuntimeError("nope")

    with pytest.raises(FallbackExhausted) as exc_info:
        client.fallback_chain("widget", [("A", always_fails), ("B", always_fails)])
    assert len(exc_info.value.attempts) == 2
    assert all(not a.success for a in exc_info.value.attempts)


def test_fallback_chain_treats_falsy_result_as_miss(tmp_path):
    client = _client(tmp_path)

    def empty():
        return None

    def real_hit():
        return {"found": True}

    result = client.fallback_chain("widget", [("empty_source", empty), ("real_source", real_hit)])
    assert result.source == "real_source"

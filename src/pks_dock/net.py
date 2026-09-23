"""
pks_dock.net
============

Single, shared HTTP layer for every external database call the pipeline
makes (NCBI, UniProt, RCSB PDB, AlphaFold DB, PubChem, ChEBI, CIR, ...).

Why this exists
----------------
Before this module, each phase script (`01_fetch_genome.py`,
`04_get_ligands.py`, `06_get_receptors.py`, `06b_predict_receptors.py`, ...)
implemented its own ad-hoc `requests.get(..., timeout=N)` calls with
inconsistent retry logic (some had none). A single flaky database should
never crash a multi-hour docking run. This module gives every phase:

1. Automatic retries with exponential backoff + jitter on transient
   failures (connection errors, timeouts, 429, 5xx).
2. On-disk response caching (`cache/http/`), keyed by URL + params, so a
   second run (or a re-run after a crash) does not re-hit rate-limited
   APIs for data that cannot have changed (e.g. a PDB entry's resolution).
3. A `fallback_chain()` helper that tries a list of (name, callable) pairs
   in order and returns the first success, logging exactly what was tried
   and why each attempt failed -- this is what lets Phase 4 walk
   PubChem -> ChEBI -> CIR -> "report compound unavailable" honestly.
4. Structured, timestamped log lines so every network decision is visible
   in the terminal and can be captured into decision_log.json by
   `pks_dock.reproducibility`.

This module does not silently swallow failures: if every source in a
fallback chain fails, the caller gets a clear `FallbackExhausted`
exception with the full attempt history, which every phase script is
expected to catch and report (never a bare crash).
"""

from __future__ import annotations

import hashlib
import json
import random
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

import requests

DEFAULT_CACHE_DIR = Path("cache/http")
DEFAULT_TIMEOUT = 30
DEFAULT_RETRIES = 3
DEFAULT_BACKOFF_BASE = 1.5  # seconds
RETRIABLE_STATUS_CODES = {429, 500, 502, 503, 504}


def _log(msg: str, phase: str = "NET") -> None:
    print(f"[{phase}] {msg}", flush=True)


@dataclass
class Attempt:
    """Record of a single fallback-chain attempt, for the decision log."""

    source: str
    success: bool
    detail: str
    elapsed_s: float


@dataclass
class FallbackResult:
    """What a successful (or exhausted) fallback_chain() call returns."""

    source: str
    value: Any
    attempts: list[Attempt] = field(default_factory=list)


class FallbackExhausted(RuntimeError):
    """Raised when every source in a fallback chain has failed."""

    def __init__(self, subject: str, attempts: list[Attempt]):
        self.subject = subject
        self.attempts = attempts
        summary = "; ".join(f"{a.source}: {a.detail}" for a in attempts)
        super().__init__(f"All sources exhausted for '{subject}'. Attempts -> {summary}")


class HTTPClient:
    """
    Cached, retrying HTTP client. One instance is safe to reuse across a
    whole phase script.

    Example
    -------
        client = HTTPClient(phase="PHASE 4")
        r = client.get_json("https://pubchem.ncbi.nlm.nih.gov/.../cids/JSON")
    """

    def __init__(
        self,
        phase: str = "NET",
        cache_dir: Path | str = DEFAULT_CACHE_DIR,
        retries: int = DEFAULT_RETRIES,
        timeout: int = DEFAULT_TIMEOUT,
        backoff_base: float = DEFAULT_BACKOFF_BASE,
        use_cache: bool = True,
        session: Optional[requests.Session] = None,
    ):
        self.phase = phase
        self.cache_dir = Path(cache_dir)
        self.retries = retries
        self.timeout = timeout
        self.backoff_base = backoff_base
        self.use_cache = use_cache
        self.session = session or requests.Session()
        if self.use_cache:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    # -- caching ----------------------------------------------------------

    def _cache_key(self, method: str, url: str, params: dict | None) -> Path:
        raw = json.dumps({"method": method, "url": url, "params": params or {}}, sort_keys=True)
        digest = hashlib.sha256(raw.encode()).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def _cache_read(self, key: Path) -> Optional[dict]:
        if not self.use_cache or not key.exists():
            return None
        try:
            return json.loads(key.read_text())
        except (json.JSONDecodeError, OSError):
            return None

    def _cache_write(self, key: Path, payload: dict) -> None:
        if not self.use_cache:
            return
        try:
            key.write_text(json.dumps(payload))
        except OSError:
            pass  # cache is a convenience, never fatal

    # -- core request with retry ------------------------------------------

    def request(
        self,
        method: str,
        url: str,
        params: dict | None = None,
        headers: dict | None = None,
        data: Any = None,
        expect_json: bool = True,
        cache: Optional[bool] = None,
    ) -> Any:
        """
        Performs a single logical request (GET/POST) with retry+backoff and
        optional disk caching. Returns parsed JSON (if expect_json) or raw
        text. Raises requests.RequestException on final failure.
        """
        use_cache = self.use_cache if cache is None else cache
        cache_key = self._cache_key(method, url, params) if use_cache and method == "GET" else None

        if cache_key is not None:
            cached = self._cache_read(cache_key)
            if cached is not None:
                _log(f"cache hit: {url}", self.phase)
                return cached["body"]

        last_exc: Exception | None = None
        for attempt in range(1, self.retries + 1):
            try:
                resp = self.session.request(
                    method, url, params=params, headers=headers, data=data, timeout=self.timeout
                )
                if resp.status_code in RETRIABLE_STATUS_CODES:
                    raise requests.HTTPError(f"HTTP {resp.status_code} (retriable)", response=resp)
                resp.raise_for_status()
                body = resp.json() if expect_json else resp.text
                if cache_key is not None:
                    self._cache_write(cache_key, {"body": body})
                return body
            except (requests.RequestException, ValueError) as exc:
                last_exc = exc
                if attempt < self.retries:
                    sleep_s = self.backoff_base * (2 ** (attempt - 1)) + random.uniform(0, 0.5)
                    _log(
                        f"attempt {attempt}/{self.retries} failed for {url} ({exc}); "
                        f"retrying in {sleep_s:.1f}s",
                        self.phase,
                    )
                    time.sleep(sleep_s)
                else:
                    _log(f"attempt {attempt}/{self.retries} failed for {url} ({exc}); giving up", self.phase)
        assert last_exc is not None
        raise last_exc

    def get_json(
        self, url: str, params: dict | None = None, headers: dict | None = None, cache: Optional[bool] = None
    ) -> Any:
        return self.request("GET", url, params=params, headers=headers, expect_json=True, cache=cache)

    def get_text(
        self, url: str, params: dict | None = None, headers: dict | None = None, cache: Optional[bool] = None
    ) -> str:
        return self.request("GET", url, params=params, headers=headers, expect_json=False, cache=cache)

    def get_binary(self, url: str, params: dict | None = None, headers: dict | None = None) -> bytes:
        """Not cached to disk as JSON (binary payloads: structures, PDFs, etc.)."""
        last_exc: Exception | None = None
        for attempt in range(1, self.retries + 1):
            try:
                resp = self.session.get(url, params=params, headers=headers, timeout=self.timeout)
                if resp.status_code in RETRIABLE_STATUS_CODES:
                    raise requests.HTTPError(f"HTTP {resp.status_code} (retriable)", response=resp)
                resp.raise_for_status()
                return resp.content
            except requests.RequestException as exc:
                last_exc = exc
                if attempt < self.retries:
                    sleep_s = self.backoff_base * (2 ** (attempt - 1)) + random.uniform(0, 0.5)
                    _log(
                        f"binary fetch attempt {attempt}/{self.retries} failed for {url} ({exc}); "
                        f"retrying in {sleep_s:.1f}s",
                        self.phase,
                    )
                    time.sleep(sleep_s)
        assert last_exc is not None
        raise last_exc

    # -- fallback chains ----------------------------------------------------

    def fallback_chain(self, subject: str, sources: list[tuple[str, Callable[[], Any]]]) -> FallbackResult:
        """
        Tries each (name, callable) in `sources` in order. `callable` takes
        no arguments, uses `self` internally (closures), and should raise
        on failure or return a falsy value to mean "no result found".

        Returns a FallbackResult on first success. Raises FallbackExhausted
        (with the full attempt history) if every source fails -- callers
        MUST catch this and report it, never suppress it silently.
        """
        attempts: list[Attempt] = []
        for name, fn in sources:
            t0 = time.monotonic()
            try:
                value = fn()
                elapsed = time.monotonic() - t0
                if value:
                    attempts.append(Attempt(name, True, "ok", elapsed))
                    _log(f"'{subject}' resolved via {name} ({elapsed:.1f}s)", self.phase)
                    return FallbackResult(source=name, value=value, attempts=attempts)
                attempts.append(Attempt(name, False, "empty result", elapsed))
                _log(f"'{subject}' via {name}: no result, falling back...", self.phase)
            except Exception as exc:  # noqa: BLE001 - deliberately broad: any source may fail any way
                elapsed = time.monotonic() - t0
                attempts.append(Attempt(name, False, str(exc), elapsed))
                _log(f"'{subject}' via {name}: FAILED ({exc}), falling back...", self.phase)
        _log(f"'{subject}': ALL sources exhausted ({len(attempts)} tried)", self.phase)
        raise FallbackExhausted(subject, attempts)

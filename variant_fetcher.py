
from __future__ import annotations

"""
Simple helper to fetch flanking genomic sequences around a variant
using the Ensembl REST API, with the UCSC Genome Browser API as a fallback.
"""

from functools import lru_cache

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ENSEMBL_REST = "https://rest.ensembl.org"
UCSC_API = "https://api.genome.ucsc.edu"
UCSC_GENOMES = {"GRCh38": "hg38", "GRCh37": "hg19"}

# (connect, read) timeouts in seconds. Ensembl can take >20 s to answer the
# first request for a region, then responds quickly once it is warmed up.
TIMEOUT = (10, 30)


def _session() -> requests.Session:
    # Keep retries short so a hanging Ensembl falls back to UCSC within ~1 min
    retry = Retry(
        total=2,
        read=1,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        respect_retry_after_header=True,
    )
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


_SESSION = _session()


def _fetch_ensembl(chrom: str, start: int, end: int, assembly: str, species: str) -> str:
    url = f"{ENSEMBL_REST}/sequence/region/{species}/{chrom}:{start}-{end}"
    params = {"coord_system_version": assembly}
    headers = {"Content-Type": "text/plain"}
    resp = _SESSION.get(url, params=params, headers=headers, timeout=TIMEOUT)
    resp.raise_for_status()
    return resp.text.strip().upper()


def _fetch_ucsc(chrom: str, start: int, end: int, assembly: str) -> str:
    genome = UCSC_GENOMES.get(assembly)
    if genome is None:
        raise ValueError(f"No UCSC genome for assembly {assembly!r}")
    # UCSC uses 0-based, half-open coordinates
    params = {"genome": genome, "chrom": f"chr{chrom}", "start": start - 1, "end": end}
    resp = _SESSION.get(f"{UCSC_API}/getData/sequence", params=params, timeout=TIMEOUT)
    resp.raise_for_status()
    return resp.json()["dna"].strip().upper()


@lru_cache(maxsize=512)
def fetch_flank_from_ensembl(
    chromosome: str,
    position: int,
    flank: int = 500,
    assembly: str = "GRCh38",
    species: str = "homo_sapiens",
) -> str:
    """
    Fetch a flanking genomic sequence around a variant using Ensembl REST.

    Falls back to the UCSC Genome Browser API (human only) if Ensembl is
    unavailable or times out. Results are cached per process.

    Args:
        chromosome: Chromosome name (e.g. "1", "2", "X")
        position: 1-based genomic position of the variant
        flank: number of bases upstream and downstream
        assembly: genome assembly, e.g. "GRCh38" or "GRCh37"
        species: species name for Ensembl, default "homo_sapiens"

    Returns:
        Flanking sequence as an uppercase string (length ~ 2*flank).
    """
    chrom = str(chromosome).replace("chr", "").replace("CHR", "")
    try:
        pos_int = int(position)
    except Exception:
        raise ValueError(f"Invalid position: {position!r}")

    start = max(1, pos_int - flank)
    end = pos_int + flank

    try:
        return _fetch_ensembl(chrom, start, end, assembly, species)
    except requests.RequestException as ensembl_err:
        if species != "homo_sapiens":
            raise
        try:
            return _fetch_ucsc(chrom, start, end, assembly)
        except (requests.RequestException, ValueError, KeyError) as ucsc_err:
            raise RuntimeError(
                "Could not fetch the flanking sequence from Ensembl or UCSC. "
                "Both services may be slow or unavailable; please try again in a minute. "
                f"(Ensembl: {ensembl_err}; UCSC: {ucsc_err})"
            ) from ensembl_err

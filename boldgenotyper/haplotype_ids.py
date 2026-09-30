"""Stable haplotype keys and display IDs.

A haplotype's stable key is seq_md5: the md5 of its uppercase, ungapped core
sequence. It is identical across runs for an identical sequence. Display IDs
(H0001, H0002, ...) carry no counts; they are numbered by member count
(descending) and then seq_md5, so identical input gives identical IDs.
"""

import hashlib
import re
from typing import Dict, Optional

_HEADER_FIELD = re.compile(r"\b(seq_md5|n_members)=(\S+)")


def seq_md5(seq: str) -> str:
    """Return the md5 hex digest of seq, uppercased with '-' and '.' removed."""
    core = str(seq).upper().replace("-", "").replace(".", "")
    return hashlib.md5(core.encode("ascii")).hexdigest()


def assign_display_ids(counts: Dict[str, int]) -> Dict[str, str]:
    """Map seq_md5 -> display ID, ordered by count descending, then seq_md5."""
    order = sorted(counts, key=lambda k: (-counts[k], k))
    width = max(4, len(str(len(order))))
    return {key: f"H{i:0{width}d}" for i, key in enumerate(order, start=1)}


def format_fasta_header(haplotype_id: str, md5: str, n_members: int) -> str:
    """Return a FASTA header line (without '>') for a haplotype."""
    return f"{haplotype_id} seq_md5={md5} n_members={int(n_members)}"


def parse_fasta_header(description: str) -> Dict[str, Optional[object]]:
    """Parse a haplotype FASTA header; fields that are absent are None."""
    text = str(description or "").lstrip(">").strip()
    parts = text.split(None, 1)
    fields = dict(_HEADER_FIELD.findall(text))
    n = fields.get("n_members")
    return {
        "haplotype_id": parts[0] if parts else None,
        "seq_md5": fields.get("seq_md5"),
        "n_members": int(n) if n is not None and n.isdigit() else None,
    }

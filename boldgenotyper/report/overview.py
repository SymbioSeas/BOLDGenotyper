"""Overview: input, key settings, counts and key findings."""

import pandas as pd

from .data import RunData
from .insights import key_findings
from .page import Section, counts, esc, kv_table, num


def _v(value, fmt="{:,}", missing="not recorded") -> str:
    return missing if value is None or (isinstance(value, float) and pd.isna(value)) else num(value, fmt)


def build(run: RunData) -> Section:
    prov, cfg, applied, s = run.provenance, run.config, run.applied, run.summary
    inp = prov.get("input") or {}
    hap, tax, coi = cfg.get("haplotype") or {}, cfg.get("taxonomy") or {}, cfg.get("coi_validation") or {}
    rooting = (applied.get("tree_rooting") or {}).get("method")
    tree = ("not built" if applied.get("n_tree_tips") is None else
            {"midpoint": "built; rooted at the midpoint", "none": "built; unrooted"}.get(
                rooting, "built; rooted on an outgroup"))
    if not applied.get("geographic_analysis"):
        regions = "no geographic analysis"
    elif applied.get("custom_shapefile"):
        regions = (f"custom shapefile <code>{esc(applied['custom_shapefile'])}</code>, field "
                   f"{esc(applied.get('shapefile_field'))} (<code>{esc(applied.get('geo_category'))}</code>)")
    else:
        regions = f"Global Oceans and Seas (<code>{esc(applied.get('geo_category'))}</code>)"
    settings = kv_table([
        ("Input", f"<code>{esc(inp.get('path'))}</code> ({_v(inp.get('n_rows'))} rows; md5 {esc(inp.get('md5'))})"),
        ("BOLD download date", esc(inp.get("download_date")) or "not recorded"),
        ("Singleton distance", _v(hap.get("min_singleton_distance"), "{}")),
        ("Species majority threshold", _v(tax.get("majority_species_threshold"), "{}")),
        ("Genetic code", _v(coi.get("mitochondrial_code"), "{}")),
        ("Tree", tree),
        ("Regions", regions),
        ("--extended", "yes" if run.extended else "no"),
    ])
    assigned = (s.get("n_assigned") or 0) + (s.get("n_assigned_corrected") or 0)
    items = [("unique records", _v(s.get("n_unique"))), ("QC-passed", _v(s.get("n_qc_pass"))),
             ("assigned", _v(assigned)), ("haplotypes", _v(s.get("n_haplotypes"))),
             ("species-level taxa", _v(s.get("n_species_haplotype_assigned")))]
    if run.region_col and not run.records.empty:
        n_regions = run.records[run.region_col].dropna().astype(str).str.strip().replace("", pd.NA).dropna().nunique()
        items.append((f"{run.region_col} values with records", _v(n_regions)))
    findings = key_findings(run)
    found = "".join(f'<li>{esc(f.text)} <a href="{esc(f.link)}">details</a></li>' for f in findings)
    html = (counts(items) + "<h3>Key findings</h3>\n"
            + (f'<ul class="findings">{found}</ul>\n' if found else '<p class="note">None.</p>\n')
            + '<p class="note">Counts only; what they mean for these taxa is for you to judge.</p>\n'
            + "<h3>Input and settings</h3>\n" + settings
            + '<p class="note">Every setting and value is in <code>run_provenance.json</code>.</p>\n')
    return Section("overview", "Overview", html)

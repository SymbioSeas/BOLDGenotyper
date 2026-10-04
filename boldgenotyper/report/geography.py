"""Geography: map overview with a species switcher, region bars, coordinate quality."""

import pandas as pd

from .data import RunData
from .page import Section, esc, kv_table, manifest_file, manifest_options, num, switcher, table


def _text(series: pd.Series) -> pd.Series:
    s = series.astype(object).where(series.notna(), "").astype(str).str.strip()
    return s.where(s != "")


def build(run: RunData) -> Section:
    rec = run.records
    html = ""
    options = manifest_options(run, "geography", "distribution_map", "All taxa")
    if options:
        links = [(t, p) for t, p in (("map (PDF)", manifest_file(run, "geography", "distribution_map", "overview_pdf")),
                                     ("one page per taxon (PDF)",
                                      manifest_file(run, "geography", "distribution_map", "pages_pdf"))) if p]
        html += ("<p>Records with geographic-quality coordinates. Choose a taxon to see its haplotypes.</p>\n"
                 + switcher("map-view", options, "Distribution map", links))
    else:
        html += '<p class="note">No map: no geographic analysis in this run, or no mapped records.</p>\n'

    for name, title in (("region_abundance_relative", "Share of each region's records"),
                        ("region_abundance_total", "Records per region")):
        opts = manifest_options(run, "extended/geography", name, "All taxa")
        if opts:
            pdf = manifest_file(run, "extended/geography", name, "pages_pdf")
            html += f"<h3>{esc(title)}</h3>\n" + switcher(f"{name}-view", opts, title,
                                                           [("one page per taxon (PDF)", pdf)] if pdf else [])

    if not rec.empty and "coord_status" in rec.columns:
        status = rec["coord_status"].astype(str)
        has_coords = rec["lat"].notna() & rec["lon"].notna() if "lat" in rec.columns else status == "ok"
        mapped = status == "ok"
        rows = [("Records", num(len(rec))), ("With coordinates", num(int(has_coords.sum()))),
                ("Mapped (geographic quality)", num(int(mapped.sum())))]
        rows += [(f"Not mapped: {k}", num(v)) for k, v in status[~mapped].value_counts().items()]
        if run.region_col:
            regions = _text(rec[run.region_col])
            rows.append((f"Mapped but outside every {run.region_col}", num(int((mapped & regions.isna()).sum()))))
        html += ("<h3>Coordinate quality</h3>\n" + kv_table(rows)
                 + '<p class="note">Country- and province-centroid, [0, 0], out-of-range and unparseable '
                   "coordinates are kept in <code>records.csv</code> but not mapped or assigned a region.</p>\n")
        if run.region_col:
            per = _text(rec[run.region_col]).dropna().value_counts().rename_axis(run.region_col).reset_index(name="records")
            html += (f"<h3>Records per {esc(run.region_col)}</h3>\n"
                     + table(per, [run.region_col, "records"], file="geography/haplotype_by_region.csv"))
    return Section("geography", "Geography", html)

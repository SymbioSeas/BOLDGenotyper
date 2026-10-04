"""Sampling and metadata: field coverage, collection dates, categorical fields."""

from .data import RunData
from .page import Section, bars, manifest_file, manifest_options, switcher, two_level_switcher

DATE_FIGURES = {"records_by_year", "sampling_timeline", "records_by_year_by_region", "metadata_coverage"}


def build(run: RunData) -> Section:
    html = "<h3>Metadata coverage</h3>\n"
    cov = run.coverage
    if not cov.empty and {"field", "pct_coverage"} <= set(cov.columns):
        items = []
        for _, r in cov.iterrows():
            pct = float(str(r["pct_coverage"]).rstrip("%") or 0)
            items.append((r["field"], pct, f"({int(r.get('n_with_value', 0)):,} of {int(r.get('n_total', 0)):,})"))
        html += bars(items, unit="%", scale=100) + '<p class="note">Share of QC-passed records with a value '
        html += "(<code>metadata_coverage.csv</code>; values normalized before counting).</p>\n"
    else:
        html += '<p class="note">No metadata coverage table.</p>\n'

    folder = "extended/metadata/figures"
    manifest = run.manifests.get(folder, {})
    if not manifest:
        html += ('<p class="note">With --extended: records by collection year, the sampling timeline '
                 "and views of each metadata field, by taxon, haplotype and region.</p>\n")
        return Section("sampling", "Sampling and metadata", html)

    year = manifest_options(run, folder, "records_by_year", "All taxa")
    by_region = manifest_options(run, folder, "records_by_year_by_region", "By region")
    if year:
        year = year[:1] + by_region + year[1:]
        html += ("<h3>Records by collection year</h3>\n"
                 '<p class="note">Counts show when specimens were sampled and deposited, not changes in '
                 "frequency.</p>\n"
                 + switcher("year-view", year, "Records by collection year",
                            [("pages (PDF)", manifest_file(run, folder, "records_by_year", "pages_pdf"))]
                            if manifest_file(run, folder, "records_by_year", "pages_pdf") else []))
    timeline = manifest_options(run, folder, "sampling_timeline", "All taxa")
    if timeline:
        html += "<h3>Sampling timeline</h3>\n" + switcher("timeline-view", timeline, "Sampling timeline")
    coverage = manifest_options(run, folder, "metadata_coverage", "All taxa")
    if coverage:
        html += "<h3>Coverage by taxon and haplotype</h3>\n" + switcher("coverage-view", coverage, "Coverage")

    fields = {}
    for name in manifest:
        if name in DATE_FIGURES or name.endswith("_by_region"):
            continue
        views = manifest_options(run, folder, name, "By taxon or haplotype")
        views[1:1] = manifest_options(run, folder, f"{name}_by_region", "By region")
        if views:
            fields[name] = views
    if fields:
        html += ("<h3>Metadata fields</h3>\n<p>Choose a field, then a view: by taxon or haplotype, by "
                 "region, or one taxon's haplotypes.</p>\n"
                 + two_level_switcher("field-view", fields, "Metadata field"))
    return Section("sampling", "Sampling and metadata", html)

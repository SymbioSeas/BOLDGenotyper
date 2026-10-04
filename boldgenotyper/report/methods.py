"""Methods: what the run did, with the values it used (run_provenance.json),
a methods statement and references."""

from pathlib import Path

import pandas as pd

from .data import RunData
from .page import Section, esc, kv_table


def _v(value, fmt="{}", missing="not recorded") -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return missing
    return fmt.format(value)


def build(run: RunData, version: str = "") -> Section:
    prov = run.provenance
    cfg, applied, tools = run.config, run.applied, prov.get("tools") or {}
    summary = applied.get("assignment_summary") or run.summary or {}
    qc, coi, core = cfg.get("qc") or {}, cfg.get("coi_validation") or {}, cfg.get("core_region") or {}
    hap, phy, tax = cfg.get("haplotype") or {}, cfg.get("phylogenetic") or {}, cfg.get("taxonomy") or {}
    notes = prov.get("data_notes") or {}
    version = prov.get("boldgenotyper_version", version)
    code = coi.get("mitochondrial_code")
    code_name = {2: "vertebrate mitochondrial", 5: "invertebrate mitochondrial"}.get(code, "NCBI genetic code")
    tree_built = applied.get("n_tree_tips") is not None
    rooting = applied.get("tree_rooting") or {}
    outgroup = esc(rooting.get("outgroup", ""))
    kept = "; support values kept on their branches"
    rooting_text = {
        "midpoint": "Midpoint: the middle of the longest path between two tips" + kept,
        "outgroup_label": f"On the outgroup tip {outgroup}" + kept,
        "outgroup_taxon": f"On the clade of haplotypes assigned to {outgroup}" + kept,
        "outgroup_fasta": f"On the outgroup sequences in {outgroup}" + kept,
        "none": "Unrooted, as written by FastTree",
    }.get(rooting.get("method"), "not recorded")
    curation = applied.get("taxon_curation") or {}
    geo = applied.get("geographic_analysis", False)
    geo_cat = applied.get("geo_category") or "ocean_basin"
    custom_shp = applied.get("custom_shapefile")
    d = hap.get("min_singleton_distance")
    threshold = tax.get("majority_species_threshold")
    inp = prov.get("input") or {}

    html = '<p class="note">What this run did, with the values it used (<code>run_provenance.json</code>).</p>\n'
    html += "<h3>Run</h3>\n" + kv_table([
        ("BOLDGenotyper", f"v{esc(version)}"),
        ("Run started", esc(_v(prov.get("started_at")))),
        ("Input", f'<code>{esc(_v(inp.get("path")))}</code> ({_v(inp.get("n_rows"), "{:,}")} rows; '
                  f'md5 {esc(_v(inp.get("md5")))})'),
        ("BOLD download date", esc(_v(inp.get("download_date")))),
        ("Command", f'<code>{esc(_v(prov.get("command_line")))}</code>'),
    ])
    html += "<h3>Quality control</h3>\n" + kv_table([
        ("Duplicate processids", f'{_v(notes.get("n_duplicate_processids"), "{:,}")} removed (first record kept)'),
        ("Orientation and reading frame",
         f"Sequences oriented and checked for an open reading frame with genetic code {_v(code)} "
         f"({code_name}): ORF coverage at least {_v(coi.get('orf_min_coverage'))}, at most "
         f"{_v(coi.get('orf_max_internal_stops'))} internal stop codons"),
        ("Length", f"At least {_v(qc.get('min_raw_length_abs'))} bp and at least "
                   f"{_v(qc.get('min_raw_length_frac_of_median'))} of the median length"),
        ("Ambiguous bases", f"N fraction at most {_v(qc.get('max_raw_N_fraction'))}"),
        ("Valid reading frame required", _v(qc.get("require_valid_orf"))),
        ("Records passing QC", f'{_v(summary.get("n_qc_pass"), "{:,}")} of '
                               f'{_v(summary.get("n_unique"), "{:,}")} unique records'),
    ])
    html += "<h3>Haplotypes</h3>\n" + kv_table([
        ("Alignment", f"QC-passed sequences aligned with MAFFT {_v(tools.get('mafft'), 'v{}', '')} (<code>--auto</code>)"),
        ("Core region", f"Alignment columns with sequence in at least {_v(core.get('core_min_coverage'))} of "
                        f"sequences; {_v(applied.get('core_region_length_bp'))} bp in this run"),
        ("Core coverage per record", f"Records whose ungapped core is shorter than {_v(core.get('min_sequence_coverage'))} "
                                     f"of the 95th-percentile core length ({_v(applied.get('min_sequence_length_bp'), '{:.0f}')} bp) "
                                     f"are excluded ({_v(summary.get('n_excluded_core'), '{:,}')} records)"),
        ("Haplotypes", "Each distinct core sequence is a haplotype (exact sequence variant), keyed by the md5 of "
                       "its sequence; haplotypes that differ only at gaps or Ns are merged"),
        ("Error filter", f"Singleton haplotypes within p-distance {_v(d)} of another haplotype are treated as "
                         "likely sequencing errors and removed"),
        ("Correction", (f"A removed singleton record is assigned to its nearest retained haplotype when that "
                        f"haplotype is unique and within {_v(d)}" if hap.get("assign_corrected", True)
                        else "Disabled (--no-assign-corrected): removed singleton records stay unassigned")),
        ("Suspect flags", f"Singletons more than {_v(hap.get('max_singleton_distance'))} from their nearest "
                          "neighbor, distance outliers and members failing the reading-frame check flag a "
                          "haplotype as suspect"),
        ("IDs", "H0001, H0002, ... in order of member count, then seq_md5"),
        ("Haplotypes retained", _v(applied.get("n_haplotypes"), "{:,}")),
    ])
    statuses = [(s, _v(summary.get(c), "{:,}")) for s, c in (
        ("assigned", "n_assigned"), ("assigned_corrected", "n_assigned_corrected"),
        ("excluded_error_singleton", "n_excluded_error"), ("excluded_core_coverage", "n_excluded_core"),
        ("failed_qc", "n_failed_qc"), ("no_sequence", "n_no_sequence"))]
    statuses.append(("Assigned (of QC-passed records)", _v(summary.get("pct_assigned"), "{}%")))
    html += "<h3>Record statuses</h3>\n" + kv_table(statuses)
    tax_rows = [
        ("Record label", "The species field, or the genus when the species is missing"),
        ("Haplotype taxon", f"The species carried by more than {_v(threshold)} of the counted labels "
                            "(members with a species-level label, plus genus-only or open labels of another "
                            "genus; members without a species label are left out; ties not resolved); "
                            "otherwise the most common genus"),
        ("Label conflicts", "Members whose species-level label differs from the haplotype’s majority label "
                            "(or from its curated taxon); listed with their evidence in label_conflicts.csv"),
        ("Divergence", "Uncorrected p-distance between aligned haplotype sequences, skipping gaps and Ns; "
                       "no taxonomic thresholds are applied"),
    ]
    if curation:
        tax_rows.append(("Curated taxa", f'{_v(curation.get("n_applied"))} haplotype(s) given the taxon in the '
                                         f'curation file <code>{esc(curation.get("path", ""))}</code> '
                                         f'(md5 {esc(_v(curation.get("md5")))}); '
                                         f'{_v(curation.get("n_not_applied"))} row(s) not applied'))
    html += "<h3>Taxonomy and labels</h3>\n" + kv_table(tax_rows)
    if tree_built:
        html += "<h3>Phylogeny</h3>\n" + kv_table([
            ("Tips", f"Haplotypes with at least {_v(applied.get('tree_min_members'))} member(s) and at least "
                     f"{_v(applied.get('tree_min_length_bp'))} bp, not suspect: {_v(applied.get('n_tree_tips'))} tips, "
                     f"{_v(applied.get('n_tree_excluded'))} excluded (phylogeny/excluded_tips.txt)"),
            ("Alignment", f"MAFFT; trimmed with trimAl ({_v(phy.get('trim_method'))}) when trimAl is installed"
             if phy.get("trim_alignment") else "MAFFT"),
            ("Inference", f"FastTree {_v(tools.get('fasttree'), 'v{}', '')} <code>-nt -gtr -gamma</code>, "
                          f"seed {_v((prov.get('seeds') or {}).get('fasttree'))}; support values are FastTree local support"),
            ("Rooting", rooting_text),
        ])
    if geo:
        ref = (f"Custom shapefile <code>{esc(Path(custom_shp).name)}</code>, field {esc(_v(applied.get('shapefile_field')))}"
               if custom_shp else "Global Oceans and Seas (GOaS) v1")
        html += "<h3>Geography</h3>\n" + kv_table([
            ("Regions", f"{ref}; region column <code>{esc(geo_cat)}</code>"),
            ("Records assigned a region", "Only records whose coordinates are of geographic quality "
                                          "(coord_status ok). Country- and province-centroid coordinates are "
                                          "kept but not assigned a region or mapped"),
            ("Coordinate status counts", ", ".join(f"{esc(k)} {n:,}" for k, n in (notes.get("coord_status") or {}).items())),
        ])
    software = [("BOLDGenotyper", f"v{esc(version)}"), ("Python", esc(_v(prov.get("python_version"))))]
    software += [(name, esc(_v(tools.get(key)))) for name, key in
                 (("MAFFT", "mafft"), ("FastTree", "fasttree"), ("trimAl", "trimal")) if tools.get(key)]
    html += "<h3>Software</h3>\n" + kv_table(software)

    n_assigned = (summary.get("n_assigned", 0) or 0) + (summary.get("n_assigned_corrected", 0) or 0)
    text = (f"COI barcode records for {esc(run.organism)} were downloaded from the Barcode of Life Data System "
            f"(BOLD; Ratnasingham &amp; Hebert, 2007) and analyzed with BOLDGenotyper v{esc(version)}. "
            f"Of {_v(summary.get('n_unique'), '{:,}', '[N]')} unique records, {_v(summary.get('n_qc_pass'), '{:,}', '[N]')} "
            f"passed quality control (length, ambiguous bases and an open reading frame under genetic code "
            f"{_v(code, missing='[N]')}). Sequences were aligned with MAFFT (Katoh &amp; Standley, 2013), and the "
            f"{_v(applied.get('core_region_length_bp'), missing='[N]')} bp core region was extracted; "
            f"{_v(summary.get('n_core_ok'), '{:,}', '[N]')} records covered it. Each distinct core sequence was "
            f"treated as a haplotype; singleton haplotypes within p-distance {_v(d, missing='[d]')} of another "
            f"haplotype were treated as likely sequencing errors, and their records assigned to the nearest "
            f"retained haplotype when unique and within that distance. {_v(n_assigned if summary else None, '{:,}', '[N]')} "
            f"records ({_v(summary.get('pct_assigned'), missing='[X]')}% of QC-passed records) were assigned to "
            f"{_v(applied.get('n_haplotypes'), '{:,}', '[M]')} haplotypes. Each haplotype was assigned the species "
            f"carried by more than {_v(threshold, missing='[t]')} of its members' counted species-level labels, "
            f"otherwise its genus")
    if curation:
        text += (f"; {_v(curation.get('n_applied'))} haplotype(s) were given a curated taxon after review of "
                 f"distance, BIN and tree evidence (curation file md5 {esc(_v(curation.get('md5')))})")
    text += ". "
    if tree_built:
        text += ("A tree of the haplotypes was inferred with FastTree (Price et al., 2010) under the GTR+Gamma model"
                 + {"midpoint": " and rooted at its midpoint. ", "none": ". "}.get(
                     rooting.get("method"), " and rooted on an outgroup. " if rooting.get("method") else ". "))
    if geo and not custom_shp:
        text += ("Records with geographic-quality coordinates were assigned to ocean basins of the Global Oceans "
                 "and Seas dataset (Flanders Marine Institute, 2021). ")
    elif geo:
        text += f"Records with geographic-quality coordinates were assigned to regions of a custom shapefile ({esc(geo_cat)}). "
    html += ("<h3>Methods statement</h3>\n"
             '<p class="note">Edit as needed; the values are this run\'s.</p>\n'
             f'<div class="statement"><p>{text}</p></div>\n')

    refs = ["Katoh, K., &amp; Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: "
            "improvements in performance and usability. <em>Molecular Biology and Evolution</em>, 30(4), 772-780.",
            "Ratnasingham, S., &amp; Hebert, P. D. (2007). BOLD: The Barcode of Life Data System. "
            "<em>Molecular Ecology Notes</em>, 7(3), 355-364."]
    if tree_built:
        refs.append("Price, M. N., Dehal, P. S., &amp; Arkin, A. P. (2010). FastTree 2 – Approximately "
                    "maximum-likelihood trees for large alignments. <em>PLoS ONE</em>, 5(3), e9490.")
        refs.append("Bianchini, G., &amp; Sánchez-Baracaldo, P. (2024). TreeViewer: Flexible, modular software to "
                    "visualise and manipulate phylogenetic trees. <em>Ecology and Evolution</em>, 14, e10873. "
                    '<a href="https://doi.org/10.1002/ece3.10873">https://doi.org/10.1002/ece3.10873</a>')
        if tools.get("trimal") and phy.get("trim_alignment"):
            refs.append("Capella-Gutiérrez, S., Silla-Martínez, J. M., &amp; Gabaldón, T. (2009). trimAl: a tool "
                        "for automated alignment trimming in large-scale phylogenetic analyses. "
                        "<em>Bioinformatics</em>, 25(15), 1972-1973.")
    if geo and not custom_shp:
        refs.append("Flanders Marine Institute (2021). Global Oceans and Seas, version 1. Available online at "
                    '<a href="https://www.marineregions.org/">https://www.marineregions.org/</a>')
    html += "<h3>References</h3>\n" + "".join(f'<p class="note">{r}</p>\n' for r in sorted(refs))
    return Section("methods", "Methods", html)

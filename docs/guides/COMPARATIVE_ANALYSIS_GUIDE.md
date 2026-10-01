# Comparing Species-Level and Family-Level Runs

This guide shows how to check BOLD species labels with BOLDGenotyper 2.0.0: first within one
run, then by comparing a run on one species with a run on its family (or genus). Every step
reads the run folders' own tables; each run's `README.md` defines their columns.

A label conflict is evidence, not a verdict. Deciding which label is right, or whether two
species really share a haplotype, is yours.

## Contents

- [What Each Comparison Shows](#what-each-comparison-shows)
- [Step 1: Run Both Analyses](#step-1-run-both-analyses)
- [Step 2: Label Conflicts Within a Run](#step-2-label-conflicts-within-a-run)
- [Step 3: Record-by-Record Crosswalk](#step-3-record-by-record-crosswalk)
- [Step 4: Records to Review](#step-4-records-to-review)
- [Interpreting the Evidence](#interpreting-the-evidence)
- [Troubleshooting](#troubleshooting)
- [Reporting the Comparison](#reporting-the-comparison)

---

## What Each Comparison Shows

**Within one run.** Every haplotype's members are listed with their labels (`labels`,
`majority_label`, `majority_fraction` in `haplotypes.csv`). A record whose species-level
label differs from its haplotype's majority label has `label_conflict` True in
`records.csv` and a row in `label_conflicts.csv`, with the evidence. A family-level run
alone finds most mislabeled records, because the records of related species are in the same
run.

**Between a species-level and a family-level run.** A species-level download holds only the
records labeled with that species, so its haplotypes cannot show that some of them share
sequences with another species. Matching the same records in a family-level run shows:

- how the species-level haplotypes map onto family-level haplotypes;
- which records of the species fall in family-level haplotypes whose majority label is
  another species.

Match records by `processid`, not haplotypes by ID or `seq_md5`. Each run aligns its own
records and extracts its own core region, so the same specimen can have a different core
sequence, haplotype ID and `seq_md5` in each run.

---

## Step 1: Run Both Analyses

Use the same `--singleton-distance` and `--genetic-code` for both runs, and a separate
output folder for each:

```bash
boldgenotyper data/Sphyrna_lewini.tsv --organism Sphyrna_lewini \
  --singleton-distance 0.015 --output Sphyrna_lewini_output

boldgenotyper data/Sphyrnidae.tsv --organism Sphyrnidae \
  --singleton-distance 0.015 --output Sphyrnidae_output
```

Choose `--singleton-distance` with `boldgenotyper-sweep` (see `PARAMETER_SWEEP_GUIDE.md`).
Download both datasets on the same day so they hold the same records.

---

## Step 2: Label Conflicts Within a Run

```python
import pandas as pd

family = "Sphyrnidae_output"
haps = pd.read_csv(f"{family}/haplotypes.csv")
conflicts = pd.read_csv(f"{family}/label_conflicts.csv")

# Haplotypes whose members carry more than one species label
mixed = haps[haps["n_species_labels"] > 1]
print(mixed[["haplotype_id", "n_members", "labels", "majority_label",
             "majority_fraction", "n_label_conflicts"]])

# One row per conflicting record, with the evidence
print(conflicts.sort_values("distance_to_nearest_with_label", ascending=False))
```

`label_conflicts.csv` columns:

| Column | Meaning |
|---|---|
| `processid`, `label`, `haplotype_id` | the record, its label and its haplotype |
| `n_members`, `n_labeled` | the haplotype's members, and those counted in the label vote |
| `majority_label`, `n_majority`, `majority_fraction` | the haplotype's most common species-level label and its share of the counted labels |
| `nearest_haplotype_with_label` | the nearest other haplotype with a member carrying the record's label (empty when none) |
| `distance_to_nearest_with_label` | p-distance from the record's haplotype to that haplotype |
| `depositor_keyword_hits`, `depositor_note` | keywords found in the depositor's notes, and the notes |

---

## Step 3: Record-by-Record Crosswalk

```python
import pandas as pd

species, family = "Sphyrna_lewini_output", "Sphyrnidae_output"
members = ["assigned", "assigned_corrected"]

sp = pd.read_csv(f"{species}/records.csv", low_memory=False)
fam = pd.read_csv(f"{family}/records.csv", low_memory=False)
fam_haps = pd.read_csv(f"{family}/haplotypes.csv").set_index("haplotype_id")

cross = (sp[sp["status"].isin(members)][["processid", "label", "haplotype_id"]]
         .merge(fam[fam["status"].isin(members)][["processid", "haplotype_id", "label_conflict"]],
                on="processid", suffixes=("_species", "_family")))
cross["family_majority"] = cross["haplotype_id_family"].map(fam_haps["majority_label"])
cross["family_majority_fraction"] = cross["haplotype_id_family"].map(fam_haps["majority_fraction"])
cross["family_labels"] = cross["haplotype_id_family"].map(fam_haps["labels"])

crosswalk = (cross.groupby(["haplotype_id_species", "haplotype_id_family",
                            "family_majority", "family_labels"], dropna=False)
             .size().rename("n_records").reset_index())
print(crosswalk.sort_values(["haplotype_id_species", "n_records"], ascending=[True, False]))
```

Each row is a pair of haplotypes (species-level, family-level) with the number of records
in both. Expected patterns:

- **One to one, or many to one.** Species-level haplotypes that differ only outside the
  family run's core region, or only at error-filtered singletons, merge into one
  family-level haplotype. Each run's core region covers the positions shared by most of
  its own records, so the two runs' core regions usually differ.
- **One to many.** A species-level haplotype split across family-level haplotypes: the
  family run's core region includes positions where its records differ.
- **Family majority is another species.** Records of the species in a family-level
  haplotype dominated by another species' label: see Step 4.

---

## Step 4: Records to Review

```python
species_label = cross["label"].mode().iloc[0]
review = cross[cross["family_majority"].notna() & (cross["family_majority"] != species_label)]
print(f"{len(review)} of {len(cross)} {species_label} records are in family-level "
      "haplotypes whose majority label is another species")

# Their evidence in the family run
conflicts = pd.read_csv(f"{family}/label_conflicts.csv")
print(conflicts[conflicts["processid"].isin(review["processid"])])
```

Records that passed QC in only one run are not in the crosswalk. To list them:

```python
only_species = set(sp.loc[sp["status"].isin(members), "processid"]) - set(cross["processid"])
print(fam[fam["processid"].isin(only_species)][["processid", "status", "status_note"]])
```

---

## Interpreting the Evidence

For each record to review, weigh:

1. **Majority strength.** `majority_fraction` near 1 with several labeled members is strong
   evidence for the majority label; a fraction near 0.5, or few labeled members, is weak.
2. **Distance to the record's own label.** A small `distance_to_nearest_with_label` means
   haplotypes carrying the record's label are close by: the species may share or nearly
   share sequences at this marker. A large distance, or none (no other haplotype carries the
   label), makes a labeling or sample error more likely.
3. **Depositor notes.** Keywords such as "cf." or "juvenile" in `depositor_keyword_hits`
   can explain an uncertain identification.
4. **Biology.** Are the species sympatric? Is hybridization or introgression known? Do the
   two species differ at this marker at all?

Records are never relabeled or removed by BOLDGenotyper. If you exclude records, list their
`processid`s and reasons in your methods or supplement.

---

## Troubleshooting

**No records in the crosswalk.** The two downloads hold different records (different dates
or search terms), or the species run's processids are not in the family download. Check
`len(set(sp.processid) & set(fam.processid))`.

**Many records missing from one run.** A record can pass QC in one run and fail core-region
coverage in the other, because each run has its own core region. Compare the `status` and
`status_note` of those records in both `records.csv` files.

**Haplotype counts differ between runs.** Expected: different core regions and different
singleton sets. Compare records, not counts.

---

## Reporting the Comparison

Suggested methods text (edit to your analysis):

> We ran BOLDGenotyper (version, citation) on the BOLD records labeled *Species* (n records)
> and on all records of *Family* (n records), downloaded on DATE, with the same singleton
> distance (value, chosen with boldgenotyper-sweep). Records were matched between runs by
> BOLD process ID. Records whose species label differed from the majority label of their
> family-level haplotype were reviewed using the majority fraction, the p-distance to the
> nearest haplotype carrying the record's label, and the depositor's notes.

Report the BOLDGenotyper version and each run's `run_provenance.json` (configuration, input
md5) so the comparison can be repeated.

---

## Related Documentation

- Tutorial 3 (`notebooks/03_comparative_analysis_quality_control.ipynb`) runs these steps
- `MIGRATION_v1_to_v2.md`: v1.0.0 output files and columns and their 2.0.0 equivalents
- Each run folder's `README.md`: every file and column

# Table 1: BOLDGenotyper Case Study Summary

| Taxon | Rank | Environment | Input records | COI sequences | COI coverage (%) | Haplotypes (full run) | Haplotypes (at elbow) | Optimal threshold | Geographic regions | Region type | Shapefile | Geo-assigned | Geo-assigned (% of COI) | Runtime (min) | Validation |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Carcharhiniformes | Order | Marine | 9244 | 6964 | 75.3 | 597 | 487 | 0.015 | 8 | ocean basins | GOaS v1 | 1145 | 16.4 | 13.5 | — |
| Sphyrnidae | Family | Marine | 1087 | 849 | 78.1 | 38 | 25 | 0.015 | 6 | ocean basins | GOaS v1 | 99 | 11.7 | 1.2 | Smith & Black (2026) |
| Panulirus | Genus | Marine | 755 | 163 | 21.6 | 87 | 53 | 0.015 | 1 | ocean basins | GOaS v1 | 1 | 0.6 | 0.8 | — |
| Salmonidae | Family | Freshwater | 5546 | 4348 | 78.4 | 276 | 260 | 0.015 | 54 | drainage basins | BasinATLAS v10 (lev. 3) | 684 | 15.7 | 5.5 | — |
| Pieridae | Family | Terrestrial | 3484 | 2345 | 67.3 | 331 | 262 | 0.085 | 7 | biogeographic realms | WWF Ecoregions 2017 | 1546 | 65.9 | 4.3 | — |


**Notes:**

- COI coverage: % of input records that contain a usable COI sequence.

- Haplotypes (full run): at default singleton-distance threshold (0.005).

- Haplotypes (at elbow): at the optimal threshold identified by parameter sweep.

- Geo-assigned: number of COI sequences spatially assigned to a geographic region.

- Geo-assigned %: proportion of COI sequences with a geographic region assignment.

- Panulirus: only 1/163 COI records carry geographic coordinates (BOLD data quality).

- Runtime on Apple M2 MacBook Pro (16 GB RAM); single-threaded pipeline.

- GOaS = General Ocean Atlas of Seas (10 named basins); BasinATLAS = HydroSHEDS basin atlas (level 3, 292 basins); WWF Ecoregions 2017 (REALM field, 8 realms).

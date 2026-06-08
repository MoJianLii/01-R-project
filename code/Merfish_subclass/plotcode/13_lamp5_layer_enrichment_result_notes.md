# Lamp5 Layer Enrichment Result Notes

Result path:

`E:/zaw/2511/mouseMerfish_zhuang_subclass/ws0.4_ss0.02/Lamp5_layer_enrichment`

Corresponding code:

`E:/zaw/01-R-project/01-R-project/code/Merfish_subclass/plotcode/13_lamp5_layer_enrichment_analysis.R`

## Overall Result

The analysis contains two complementary levels:

1. Cluster-level enrichment: based on `mouse_subclass_cluster_total_over3.txt`, which keeps merged enriched regions with `enrich_subclass_cell_ids_num >= 3`.
2. Raw cell-level distribution: based on all mapped `neocortex_new` cells, using `Merfish_mouse_neocortex_layer_region.txt` to assign each cell to cortical layer.

The main conclusion is that Lamp5 inhibitory neurons show the strongest enrichment and raw-cell representation in layer 1, with a secondary but weaker presence in layer 2/3. Layers 4, 5, and 6a have fewer Lamp5 cells and fewer Lamp5-enriched clusters. Layer 6b contains some raw Lamp5 cells but no Lamp5 enriched clusters after the over3 region filter.

## Key Numbers

### Cluster-Level Enrichment

| Layer | Lamp5 clusters | GABA clusters | Lamp5 / GABA clusters | Lamp5 enriched cells |
|---|---:|---:|---:|---:|
| 1 | 382 | 641 | 59.6% | 2,265 |
| 2/3 | 317 | 1,605 | 19.8% | 1,855 |
| 4 | 42 | 514 | 8.2% | 161 |
| 5 | 127 | 1,220 | 10.4% | 560 |
| 6a | 100 | 805 | 12.4% | 393 |
| 6b | 0 | 16 | 0% | 0 |

Cluster-level interpretation:

Layer 1 is the dominant Lamp5-enriched layer by both absolute cluster count and relative proportion among GABA clusters. Layer 2/3 has many Lamp5 clusters in absolute count, but because the total number of GABA clusters is high, its relative Lamp5 fraction is much lower than layer 1.

### Raw Cell-Level Distribution

| Layer | Raw cells | GABA cells | Lamp5 cells | Lamp5 / all cells | Lamp5 / GABA cells |
|---|---:|---:|---:|---:|---:|
| 1 | 153,268 | 8,289 | 4,536 | 2.96% | 54.7% |
| 2/3 | 319,304 | 20,669 | 4,003 | 1.25% | 19.4% |
| 4 | 103,176 | 6,810 | 542 | 0.53% | 8.0% |
| 5 | 289,609 | 26,074 | 1,455 | 0.50% | 5.6% |
| 6a | 218,194 | 12,408 | 1,026 | 0.47% | 8.3% |
| 6b | 20,135 | 782 | 151 | 0.75% | 19.3% |

Raw-cell interpretation:

Layer 1 again shows the strongest Lamp5 signal: more than half of GABA cells in layer 1 are Lamp5-classified cells. Layer 2/3 also has a large absolute number of Lamp5 cells, but the Lamp5/GABA fraction is much lower than layer 1. The 6b fraction should be interpreted cautiously because the denominator is small and no over3 Lamp5 enriched cluster survives in 6b.

## Figure-by-Figure Interpretation

### Fig_Lamp5_01_cluster_count_by_layer

This bar plot shows the number of Lamp5-enriched merged regions per layer.

Main message: Lamp5 enriched clusters are concentrated in layer 1 and layer 2/3. Layer 1 has the highest number of Lamp5 clusters, and layer 6b has none after the over3 filter.

### Fig_Lamp5_02_cluster_fraction_among_gaba_by_layer

This plot normalizes Lamp5 clusters by all GABA clusters in the same layer.

Main message: layer 1 is not only high in absolute count but also highly specific. About 59.6% of layer 1 GABA enriched clusters are Lamp5 clusters, much higher than layer 2/3 at about 19.8%.

### Fig_Lamp5_03_raw_cell_count_by_layer

This plot shows raw Lamp5 cell counts from the original neocortex files.

Main message: raw Lamp5 cells are most abundant in layer 1 and layer 2/3. This supports that the cluster-level enrichment is not an artifact of region merging alone.

### Fig_Lamp5_04_raw_cell_fraction_among_gaba_by_layer

This plot shows Lamp5 cells as a fraction of all GABA cells per layer.

Main message: layer 1 is the strongest Lamp5-dominant inhibitory layer. Layer 2/3 and 6b have similar fractions around 19%, but 6b has far fewer total GABA cells and no surviving over3 Lamp5 cluster, so it is less robust.

### Fig_Lamp5_05_cluster_size_vs_other_gaba

This plot compares Lamp5-enriched cluster sizes against other GABA clusters by layer.

Main message: Lamp5 clusters are not just counted; their merged-region sizes can be compared with the broader GABA cluster background. This is useful to check whether Lamp5 enrichment is driven by unusually large or small merged regions.

### Fig_Lamp5_06_lamp5_region_ei_ratio_by_layer

This plot shows E/I ratio inside Lamp5-enriched merged regions:

`E/I ratio = GLU cell number / GABA cell number`

Main message: Lamp5-enriched regions are embedded in different excitatory/inhibitory environments across layers. Layer 1 has median E/I ratio near 1, while deeper layers and layer 2/3 have higher median E/I ratios, suggesting Lamp5-enriched regions in those layers are surrounded by more excitatory cells relative to GABA cells.

### Fig_Lamp5_07_subclass_composition_by_layer

This plot separates Lamp5-labelled subclasses.

Main message: almost all Lamp5 signal is `049 Lamp5 Gaba`. A small number of `050 Lamp5 Lhx6 Gaba` clusters appears in layer 5 and layer 6a only.

### Fig_Lamp5_08_raw_spatial_top_slices

This plot shows representative slices with the largest raw Lamp5 cell counts.

Main message: Lamp5 cells are spatially visible as enriched bands or local fields rather than being uniformly scattered. Most top raw Lamp5 slices come from layer 1 or layer 2/3, consistent with the layer-level summaries.

### Fig_Lamp5_09_slice_layer_variability

This plot shows slice-to-slice variation in Lamp5/GABA fraction by layer.

Main message: layer 1 has a high and broad Lamp5 fraction across many slices, while deeper layers show lower fractions and more dependence on individual slice-layer composition.

## Important Caveats

Layer 6b should not be overinterpreted. Raw cell counts detect 151 Lamp5 cells in 6b, but only 782 total GABA cells are present across mapped 6b, and no Lamp5 enriched cluster remains after the over3 cluster filter.

Cluster-level and raw-cell-level results answer different questions:

Cluster-level results indicate where Lamp5 forms significant merged enriched regions.

Raw-cell-level results indicate where Lamp5 cells exist anatomically, regardless of whether they form a merged enriched region.

Therefore, the strongest and most reliable conclusion is the concordant layer 1 dominance, because it is supported by both cluster-level enrichment and raw-cell distribution.

## Suggested Manuscript-Style Summary

Lamp5 inhibitory neurons displayed a strong layer-specific spatial signature in the neocortex. At the merged enriched-region level, Lamp5 clusters were most frequent in layer 1 and accounted for nearly 60% of all GABA-enriched clusters in this layer. Raw-cell analysis confirmed this pattern, with Lamp5 cells comprising 54.7% of layer 1 GABA cells, compared with 19.4% in layer 2/3 and below 10% in layers 4, 5, and 6a. These results indicate that Lamp5 interneurons are preferentially organized in layer 1-enriched spatial domains, while their deeper-layer occurrences are less frequent and less consistently organized into enriched merged regions.

# Raw Spatial Data

10X Visium spatial transcriptomics data for mouse intervertebral disc (IVD) mutant vs control comparison on slide **V11A19-364**, capture area **D1**, library **IVD2**.

All samples below use the same sequencing FASTQs (`IVD2_S2_*`); they differ only in manual Loupe tissue alignment (ROI masks).

## IVD_only

Strict IVD-only tissue masks (`Load10X_Spatial`).

| Role | Space Ranger output | Loupe alignment | Spots |
|------|---------------------|-----------------|-------|
| Mutant | `spaceranger_outs/Mutant_IVD_manual/outs/` | `loupe_alignment/Mutant_IVD_V11A19-364-D1.json` | 360 |
| Control | `spaceranger_outs/Control_IVD_manual/outs/` | `loupe_alignment/Control_IVD_V11A19-364-D1.json` | 371 |

Tissue masks exclude bone and non-IVD regions.

## full_alignment_with_bone

Broader manual alignments including bone-adjacent tissue.

| Role | Space Ranger output | Loupe alignment | Spots |
|------|---------------------|-----------------|-------|
| Mutant | `spaceranger_outs/M_IVD2_P1_manual/outs/` | `loupe_alignment/Mutant_V11A19-364-D1.json` | 827 |
| Control | `spaceranger_outs/C1_IVD2_P1_manual/outs/` | `loupe_alignment/Control_wo_Neuronal_V11A19-364-D1.json` | 1063 |

Control mask excludes neuronal regions but retains bone-associated spots.

## Per-sample `outs/` contents

Each `outs/` folder contains:

- `filtered_feature_bc_matrix.h5`
- `raw_feature_bc_matrix.h5`
- `raw_probe_bc_matrix.h5` (where available)
- `spatial/` (`tissue_positions.csv`, `scalefactors_json.json`, tissue images)
- `metrics_summary.csv`
- `probe_set.csv`

## Loading in Seurat (R)

```r
library(Seurat)

# IVD-only ROI
mutant  <- Load10X_Spatial("raw_data/IVD_only/spaceranger_outs/Mutant_IVD_manual/outs")
control <- Load10X_Spatial("raw_data/IVD_only/spaceranger_outs/Control_IVD_manual/outs")

# Full alignment with bone
mutant_full  <- Load10X_Spatial("raw_data/full_alignment_with_bone/spaceranger_outs/M_IVD2_P1_manual/outs")
control_full <- Load10X_Spatial("raw_data/full_alignment_with_bone/spaceranger_outs/C1_IVD2_P1_manual/outs")
```

## Source

Processed on the Gray lab server from Space Ranger count outputs under `/stor/work/Gray/2022Spatial_02.count/`.

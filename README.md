# HiRepliMap

Annotation of high-resolution Repli-seq features: Applicable to Repli-seq data with more than 4 S-phase fractions.

# Overview
This script is designed to annotate features of High-resolution Repli-seq.

# Features Identified
Below are the features identified in the analysis:

- **Initiation Zone (IZs):** The initiation of replication sites. It is identified as the binwise max repli-seq values flanked by consecutive bins with max signals in the later S fractions.
- **Right or Left TTR (R_TTR / L_TTR):** The transition time region is identified as consecutive bins with gradually increased or decreased S fractions with max signals, robust to the regions with max signals staying in the same regions smaller than **n** bins.
- **Termination Zone (TZs):** The termination of replication sites. It is identified as the binwise max repli-seq values flanked by consecutive bins with max signals in the earlier S fractions.
- **Late Constant Time Region (Late CTRs):** The CTR region is identified as the consecutive bins >=4 which the max signals stay in the same S fractions (restricted to S14-S16 fractions) and flanked by consecutive bins with max signals in the earlier S fractions.
- **Steady breaks (SBs):** The max signals stay in the same fractions and do not belong to other features.
- **Non determined (ND):** Bins without signals.

# Data Processing
Please follow the **Preprocess workflow** to obtain `repliseq_normalized.igv` files.

# Running the Script
## Input Parameters:
```bash
-i : repliseq_normalized.igv   # IGV files from the preprocess workflow
-o : the path for the output folder # output will be RepliFeatures.csv 
-n : the number of bins that remain within the same S-phase fractions and are tolerated as TTRs 
-s : how many s frations in repliseq_normalized.igv 
```

# Final output

- RepliFeatures.csv #The final, corrected annotation is available in the manually_annot2 column, or you can simply use the *_merge.bed files.
- *_merge.bed 

# Example

```bash

module load BEDTools/2.31.1-GCC-14.1.0 # The script requires Bedtools
module load Python/3.12.4-GCCcore-14.1.0

python 1_Annot_features.py -i /omics/groups/OE0574/internal/boyu/result_repliseq_lichin/repliseq_normalized.igv -o /home/l538g/workingf/brainbreaks/DSB/Repliseq_smooth/github/Annotation-of-high-resolution-Repli-seq-features/Test -n 5 -s 16
```

# For Script QC
This compares annotated replication features between two Repli-seq replicates and generates QC visualizations (or two conditions as you like).


## Inputs:

- rep1, rep2: CSVs with annotated Repli-seq features (two files you want to compare).

- rep1_name, rep2_name: Labels for the replicates or conditions.

- sample: Sample name used in output files.

## out: Output directory.

Outputs:
*_Percentage_compare.pdf: Stacked bar plot showing the feature composition per replicate.

*_Percentage_shift.pdf: Heatmap showing feature annotation consistency between replicates. 
(It compares the distribution of features between repeat1 and repeat2. Each cell represents the number of elements that are annotated as a given feature in repeat1 and simultaneously as another (or the same) feature in repeat2. For example, among all the elements labeled as IZ in repeat2, what proportion were IZ, CTR, or other feature types in repeat1. Each column sum is 1. It captures how features are retained or reclassified between repeat1 and repeat2.)

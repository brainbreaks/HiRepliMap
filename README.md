# Annotation-of-high-resolution-Repli-seq-features
16 fractions of Repli-seq data


# Overview
This script is designed to annotate features of High-resolution Repli-seq.

# Features Identified
Below are the features identified in the analysis:

- **Initiation Zone (IZs):** The initiation of replication sites. It is identified as the binwise max repli-seq values flanked by consecutive bins with max signals in the later S fractions.
- **Right or Left TTR (R_TTR / L_TTR):** The transition time region is identified as consecutive bins with gradually increased or decreased S fractions with max signals, robust to the regions with max signals staying in the same regions smaller than 3 bins.
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
-o : the path for the output file and name # output will be ~path/*name_rep1RepliFeatures.csv 

```

# Example

```bash
python 1_Annot_features.py -i /omics/groups/OE0574/internal/boyu/result_repliseq_lichin/repliseq_normalized.igv -o /home/l538g/workingf/brainbreaks/DSB/Repliseq_smooth/github/Annotation-of-high-resolution-Repli-seq-features/Test
```


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

python 1_Annot_features_v3.py -i /omics/groups/OE0574/internal/boyu/result_repliseq_lichin/repliseq_normalized.igv -o /home/l538g/workingf/brainbreaks/DSB/Repliseq_smooth/github/Annotation-of-high-resolution-Repli-seq-features/Test -n 5 -s 16
```



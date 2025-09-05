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


*_Percentage_compare.pdf: Stacked bar plot showing the feature composition per replicate.

*_Percentage_shift.pdf: Heatmap showing feature annotation consistency between replicates. 
(It compares the distribution of features between repeat1 and repeat2. Each cell represents the number of elements that are annotated as a given feature in repeat1 and simultaneously as another (or the same) feature in repeat2. For example, among all the elements labeled as IZ in repeat2, what proportion were IZ, CTR, or other feature types in repeat1. Each column sum is 1. It captures how features are retained or reclassified between repeat1 and repeat2.)

------------------------------------------------------------------------------
# Preprocess workflow

The data processing follows the methods outlined in this paper:
Zhao, P.A., Sasaki, T. & Gilbert, D.M. High-resolution Repli-Seq defines the temporal choreography of initiation, elongation and termination of replication in mammalian cells. Genome Biol 21, 76 (2020). https://doi.org/10.1186/s13059-020-01983-8

1. Reads per million (RPM) were computed using a 50-kb bin size
2. Each S phase fraction was normalized to G1 using the log2 ratio
3. Data underwent smoothing through the application of a Gaussian filter
4. Each bin was scaled to 100

### For new experiments **without** G1_bin50000.bdg [**Start here**]

Modify the file path according to your own data

#### Input Parameters
- -c : Specify the path of the file with chromosome information [Chromosome size such as ChromInfo-hg19.txt ]
- -b: Provide the path to the G1 bam file
- -o: Set the path for your output files

```bash
bsub -R "rusage[mem=200G]" -q long -n 40 bash Produce_G1_table.sh -c ~/ChromInfo-hg19.txt -b ~/DMSO_G1.bam -o ~/test
```

#### Outputs


1. G1_bin50000.bdg  ## **The g1 input for the repliseq_normalization.py**
2. G1_ref_bin50000.bed  
3. G1_ref_bin50000_sorted.bed  
4. G1_sorted.bam


### For new experiments **with** G1 data [**Start here**]

```bash
python repliseq_normalization.py --help
```

#### Input Parameters


- -i : Specify the path of the input folder containing your Repli-seq data [S1 to S16 bin50000.bdg files]
- -g1: Provide the path to the G1_bin50000.bdg file
- -o: Set the path for your output files

```
- S1 to S16 bin50000.bdg files generated by Sergej's script
- Ensure that the bin50000.bdg files for S1 to S16 follow the naming pattern S*_bin50000.bdg

```

```bash
python repliseq_normalization.py -i ~/bedgraph -g1 ~/G1_bin50000.bdg -o Result

```

#### Outputs

1. G1_repliseq.csv  # G1 RPM data
2. repliseq_Sfractions.csv # Merge S fraction data
3. repliseq_normalized.igv
4. S1 to S16 bedGraph files

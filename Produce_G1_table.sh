#!usr/bin/bash


echo "load tools"
module load bedtools/2.29.2
module load samtools/1.9


while getopts c:b:o: flag
do
    case "${flag}" in
        c) chromoinfo=${OPTARG};;
        b) g1_bam=${OPTARG};;
        o) output=${OPTARG};;
    esac
done
echo "chromosome information: $chromoinfo";
echo "G1 BAM file: $g1_bam";
echo "Output folder: $output";

bedtools makewindows -w 50000 -g $chromoinfo > $output/G1_ref_bin50000.bed

sort -k1,1 -k2,2n $output/G1_ref_bin50000.bed > $output/G1_ref_bin50000_sorted.bed

echo "Sort G1 BAM file ... may take some time"

samtools sort -@ 40 -o $output/G1_sorted.bam $g1_bam

echo "Generate G1_bin50000 bedgraph"
bedtools intersect -abam $output/G1_ref_bin50000_sorted.bed -b $output/G1_sorted.bam -c -bed -g $chromoinfo -sorted > $output/G1_bin50000.bdg

echo "Done"

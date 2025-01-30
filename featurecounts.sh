#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --cpus-per-task=8
#SBATCH --mem=10000MB
#SBATCH --time=14:00:00
#SBATCH --output=/home/xzhou/rnaseq_course/counts/featurecounts_%j.out
#SBATCH --error=/home/xzhou/rnaseq_course/counts/featurecounts_%j.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xin.zhou@students.unibe.ch

module load apptainer

MAPPED_DIR="/home/xzhou/rnaseq_course/mapping_result"
COUNT_DIR="/data/users/xzhou/rnaseq_course/counts"
ANNOTATION="/home/xzhou/reference_gene/Homo_sapiens.GRCh38.113.gtf"
mkdir -p $COUNT_DIR

# FeatureCounts for individual samples
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_HER21.txt -p $MAPPED_DIR/HER21_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_HER22.txt -p $MAPPED_DIR/HER22_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_HER23.txt -p $MAPPED_DIR/HER23_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_NonTNBC1.txt -p $MAPPED_DIR/NonTNBC1_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_NonTNBC2.txt -p $MAPPED_DIR/NonTNBC2_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_NonTNBC3.txt -p $MAPPED_DIR/NonTNBC3_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_Normal1.txt -p $MAPPED_DIR/Normal1_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_Normal2.txt -p $MAPPED_DIR/Normal2_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_Normal3.txt -p $MAPPED_DIR/Normal3_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_TNBC1.txt -p $MAPPED_DIR/TNBC1_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_TNBC2.txt -p $MAPPED_DIR/TNBC2_mapped.bam
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_TNBC3.txt -p $MAPPED_DIR/TNBC3_mapped.bam
# FeatureCounts for all samples together (for DESeq2)
featureCounts -a $ANNOTATION -o $COUNT_DIR/output_counts_all_samples.txt -p \
  $MAPPED_DIR/HER21_mapped.bam $MAPPED_DIR/HER22_mapped.bam $MAPPED_DIR/HER23_mapped.bam \
  $MAPPED_DIR/NonTNBC1_mapped.bam $MAPPED_DIR/NonTNBC2_mapped.bam $MAPPED_DIR/NonTNBC3_mapped.bam \
  $MAPPED_DIR/Normal1_mapped.bam $MAPPED_DIR/Normal2_mapped.bam $MAPPED_DIR/Normal3_mapped.bam \
  $MAPPED_DIR/TNBC1_mapped.bam $MAPPED_DIR/TNBC2_mapped.bam $MAPPED_DIR/TNBC3_mapped.bam

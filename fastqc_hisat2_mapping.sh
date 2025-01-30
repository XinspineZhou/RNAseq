#!/bin/bash

#SBATCH --job-name=rnaseq_pipeline
#SBATCH --cpus-per-task=10
#SBATCH --mem=8000MB
#SBATCH --time=20:00:00
#SBATCH --output=/home/xzhou/rnaseq_course/logs/rnaseq_pipeline_%j.out
#SBATCH --error=/home/xzhou/rnaseq_course/logs/rnaseq_pipeline_%j.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xin.zhou@students.unibe.ch

# Define directories
READS_DIR="/data/courses/rnaseq_course/breastcancer_de/reads"
OUTPUT_DIR="/home/xzhou/rnaseq_course"
FASTQC_DIR="$OUTPUT_DIR/fastqc_results"
MAPPING_DIR="$OUTPUT_DIR/mapping_result"
GENOME_DIR="/scratch/xzhou/reference_gene"
GENOME_FASTA="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOME_GTF="$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf"

# Create necessary directories
mkdir -p $FASTQC_DIR $MAPPING_DIR/sam_to_bam $MAPPING_DIR/sort_bam $GENOME_DIR $OUTPUT_DIR/logs

# Step 1: Run FastQC on all samples
echo "Step 1: Running FastQC..."
module load FastQC/0.11.9-Java-11
for file in $READS_DIR/*_R1.fastq.gz; do
    sample=$(basename $file _R1.fastq.gz)
    fastqc -o $FASTQC_DIR $READS_DIR/${sample}_R1.fastq.gz $READS_DIR/${sample}_R2.fastq.gz
done

# Step 2: Download Reference Genome and Annotation
echo "Step 2: Downloading reference genome..."
wget -P $GENOME_DIR https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -P $GENOME_DIR https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip $GENOME_DIR/*.gz

# Step 3: Run Hisat2 Mapping, SAM to BAM, Sorting, and Indexing
echo "Step 3: Running Hisat2 mapping..."
samples=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")

for sample in "${samples[@]}"; do
    echo "Processing $sample..."
    
    # Run Hisat2 mapping
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    hisat2 -p 4 -x $GENOME_DIR/GRCh38_index \
        -1 $READS_DIR/${sample}_R1.fastq.gz \
        -2 $READS_DIR/${sample}_R2.fastq.gz \
        -S $MAPPING_DIR/${sample}_mapped.sam"

    # Convert SAM to BAM
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    samtools view -@ 1 -bS $MAPPING_DIR/${sample}_mapped.sam > $MAPPING_DIR/${sample}_mapped.bam"
    rm $MAPPING_DIR/${sample}_mapped.sam

    # Sort BAM file
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    samtools sort -@ 4 -m 2000M \
        -T $MAPPING_DIR/temp_sort_${sample} \
        -o $MAPPING_DIR/${sample}_sorted.bam \
        $MAPPING_DIR/${sample}_mapped.bam"

    # Index BAM file
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    samtools index $MAPPING_DIR/${sample}_sorted.bam"

    echo "Completed processing for sample $sample."
done

echo "RNA-Seq processing pipeline completed!"

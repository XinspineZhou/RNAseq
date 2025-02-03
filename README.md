# RNA-seq Analysis Pipeline

This pipeline processes RNA sequencing data, including mapping, differential expression analysis, and gene ontology (GO) enrichment. The pipeline consists of four main steps:

1. **FastQC, Read Mapping with Hisat2 andGene Quantification with FeatureCounts**
2. **Gene Ontology Enrichment (GOE)**
3. **Differential Expression Analysis with DESeq2**
4. **GO Enrichment Analysis with GoEnrich**

---
## 1. FastQC, Read Mapping with Hisat2 andGene Quantification with FeatureCounts
Step 1: Quality Control with FastQC

Before proceeding with mapping, FastQC is used to assess the quality of raw sequencing reads.

Step 2: Read Mapping with HISAT2
    •    Input: Paired-end FASTQ files (*_R1.fastq.gz, *_R2.fastq.gz).
    •    Reference Genome: GRCh38 (Indexed with HISAT2).
    •    Output: Sorted and indexed BAM files.

Step 3: Gene Quantification with FeatureCounts
    •    Input: HISAT2-aligned BAM files.
    •    Annotation File: Homo_sapiens.GRCh38.113.gtf.
    •    Output:
    •    Individual count files for each sample.
    •    A combined count matrix (output_counts_all_samples.txt) for DESeq2 analysis.
    
Usage:  sbatch fastqc_hisat2_mapping.sh
        sbatch featurecounts.sh
        
Script Highlights:
    •    Uses apptainer to run HISAT2 within a Singularity container.
    •    Converts SAM to BAM, sorts, and indexes BAM files using samtools.
    •    Uses FeatureCounts to assign reads to annotated genes.
    •    Generates gene count tables for differential expression analysis.

---

## 2. Differential Expression Analysis with DESeq2

The `DESeq2.R` script conducts differential expression analysis.

### **Usage:**
Run the script in R:

```r
source("DESeq2.R")
```

### **Key Features:**
- Uses `DESeq2` to identify differentially expressed genes.
- Compares **HER2, NonTNBC, TNBC vs. Normal** conditions.
- Outputs tables of **significant genes**.

---

## 3. GO Enrichment Analysis with GoEnrich

The `GoEnrich.R` script specifically handles GO enrichment for differentially expressed genes.

### **Usage:**
Run the script in R:

```r
source("GoEnrich.R")
```

### **Key Features:**
- **Extracts** significant genes from DESeq2.
- **Performs GO enrichment analysis** using `enrichGO()`.
- **Removes redundant calculations** found in `GOE.R` for efficiency.
- Outputs **GO enrichment tables** and **visualizations**.

---

## Notes:
- **Ensure all required R packages are installed** (`ggplot2`, `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, etc.).
- **Check for redundant code** in `GoEnrich.R` and `GOE.R`. GOE analysis should be consolidated into `GoEnrich.R`.
- **Results** will be stored in output directories created during execution.

---

## Citation:
If using this pipeline, cite the following:
- **HISAT2**: Kim et al., 2015
- **DESeq2**: Love et al., 2014
- **clusterProfiler**: Yu et al., 2012

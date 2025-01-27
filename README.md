# RNAseq
#!/bin/bash
#SBATCH --job-name=hisat2_mapping
#SBATCH --cpus-per-task=10             # 使用10个CPU核心
#SBATCH --mem=8000MB                  # 每个任务分配8GB内存
#SBATCH --time=20:00:00               # 最长运行时间为13小时
#SBATCH --output=/home/xzhou/rnaseq_course/mapping_result/hisat2_mapping_%j.out
#SBATCH --error=/home/xzhou/rnaseq_course/mapping_result/hisat2_mapping_%j.err
#SBATCH --partition=pibu_el8          # 使用指定的分区
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xin.zhou@students.unibe.ch

# 定义参考基因目录和样本列表
REF_DIR="/home/xzhou/reference_gene"
samples=("HER22" "HER23" "NonTNBC1" "Normal1" "TNBC1" "NonTNBC2" "Normal2" "TNBC2" "NonTNBC3" "Normal3" "TNBC3")

# 创建输出目录
mkdir -p /home/xzhou/rnaseq_course/mapping_result/sam_to_bam
mkdir -p /home/xzhou/rnaseq_course/mapping_result/sort_bam

# 循环处理每个样本
for sample in "${samples[@]}"
do
    echo "Processing sample $sample..."

    # Step 1: Hisat2 mapping
    echo "Running Hisat2 for sample $sample..."
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    hisat2 -p 4 -x $REF_DIR/GRCh38_index \
        -1 /data/courses/rnaseq_course/breastcancer_de/reads/${sample}_R1.fastq.gz \
        -2 /data/courses/rnaseq_course/breastcancer_de/reads/${sample}_R2.fastq.gz \
        -S /home/xzhou/rnaseq_course/mapping_result/${sample}_mapped.sam"

    # Step 2: Convert SAM to BAM
    echo "Converting SAM to BAM for sample $sample..."
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    samtools view -@ 1 -bS /home/xzhou/rnaseq_course/mapping_result/${sample}_mapped.sam > \
        /home/xzhou/rnaseq_course/mapping_result/${sample}_mapped.bam"
    rm /home/xzhou/rnaseq_course/mapping_result/${sample}_mapped.sam

    # Step 3: Sort BAM file
    echo "Sorting BAM file for sample $sample..."
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    samtools sort -@ 4 -m 2000M \
        -T /home/xzhou/rnaseq_course/mapping_result/temp_sort_${sample} \
        -o /home/xzhou/rnaseq_course/mapping_result/${sample}_sorted.bam \
        /home/xzhou/rnaseq_course/mapping_result/${sample}_mapped.bam"

    # Step 4: Index BAM file
    echo "Indexing BAM file for sample $sample..."
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    samtools index /home/xzhou/rnaseq_course/mapping_result/${sample}_sorted.bam"

    echo "Completed processing for sample $sample."
done

echo "All samples processed sequentially."

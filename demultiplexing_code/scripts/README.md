# Amplicon Sequencing Pipeline — Galaxy-compatible

## Tools

| Step | CLI tools | Galaxy tools | Input | Output |
|------|---------|-----------|------|------|
| 1 | `fastqc` | FastQC | raw FASTQ | QC HTML |
| 2 | `fastp` | fastp | R1+R2 FASTQ | merged FASTQ |
| 3 | `demux_by_index.py` | Je Demultiplex | merged + barcode file | per sample FASTQ |
| 4 | `cutadapt --demultiplex` | Cutadapt | per sample FASTQ + barcodes.fasta | per amplicon FASTQ |
| 5 | `bwa mem` | Map with BWA-MEM | amplicon FASTQ + ref | BAM |
| 6 | `samtools idxstats` | SAMtools idxstats | BAM | reads count matrix |
| 7 | `multiqc` | MultiQC | all QC output | HTML report |

---

## Quick start

```bash
# 1. install dependencies
conda create -p /hpc/group/youlab/xc232/reads_demultiplexed -c bioconda -c conda-forge fastqc fastp cutadapt bwa samtools multiqc pandas pysam -y

# To activate this environment, use
#
#     $ conda activate /hpc/group/youlab/xc232/reads_demultiplexed
#
# To deactivate an active environment, use
#
#     $ conda deactivate

# 2. construct amplicon references

$ conda activate /hpc/group/youlab/xc232/reads_demultiplexed

$ cd reads_05182026

python3 scripts/build_amplicon_reference.py \
    --primers-csv sink_dual_amplicons.csv \
    --out-fasta reference/amplicons.fa

# 3. pipeline_config.sh
chmod +x scripts/amplicon_pipeline.sh
. amplicon_pipeline.sh all
```

---

### Index files（Plate1_SampleSheet.txt）
```
Sample_ID   I7_Index_ID   index      I5_Index_ID   index2
1-1         D701          CGAGTAAT   D501          TATAGCCT
```
Step 3会自动转换为Je barcode格式：`1-1:CGAGTAAT:TATAGCCT`

### Barcode + Amplicon文件（primers CSV）
```
Seq_ID, Target_Barcode, Forward_Primer, Reverse_Primer, Amplicon_Sequence, ...
```
- Step 4用`Target_Barcode`+`Forward_Primer`构建cutadapt linked adapter FASTA
- Step 5用`Amplicon_Sequence`构建比对参考

### 在Galaxy中使用自定义输入
1. 上传`Plate1_SampleSheet.txt` → 运行`samplesheet_to_je_barcodes.py`生成barcode文件 → 上传给Je Demultiplex
2. 上传primers CSV → 运行`primers_to_cutadapt_barcodes.py`生成barcodes.fasta → 上传给Cutadapt
3. 上传primers CSV → 运行`build_amplicon_reference.py`生成amplicons.fa → 上传给BWA-MEM

---

## 输出文件

| 文件 | 说明 |
|------|------|
| `results/06_counts/count_matrix.tsv` | **主输出**：样本×扩增子reads计数矩阵 |
| `results/06_counts/count_matrix_long.tsv` | 长格式，适合R/Python导入 |
| `results/06_counts/amplicon_summary.tsv` | 每个扩增子的统计摘要 |
| `results/06_counts/sample_totals.tsv` | 每个样本总mapped reads |
| `results/07_report/multiqc_report.html` | 综合QC报告 |

### count_matrix.tsv格式
```
Seq_ID          1-1    1-2    1-3    ...
A. baylyi ADP1  1523   1847   902    ...
DA28102_1       3201   2984   3456   ...
...
```

---

## Pipeline流程图

```
原始FASTQ (R1+R2)
    │
    ▼ Step 1: FastQC (质控)
    │
    ▼ Step 2: fastp (adapter修剪 + PE合并 → merged.fastq.gz)
    │
    ▼ Step 3: Index demux (SampleSheet → 96个样本FASTQ)
    │         Je Demultiplex / demux_by_index.py
    │
    ▼ Step 4: Barcode demux (Target_Barcode → 每样本×89扩增子)
    │         Cutadapt --demultiplex / linked adapter
    │
    ▼ Step 5: BWA-MEM (比对到amplicons.fa参考)
    │
    ▼ Step 6: SAMtools idxstats (每个BAM → reads计数)
    │         merge_idxstats.py → count_matrix.tsv
    │
    ▼ Step 7: MultiQC (汇总报告)
```

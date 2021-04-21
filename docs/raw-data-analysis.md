# Detecting SARS-CoV-2 lineages and mutational load in municipal wastewater; raw data analysis


## Steps


### Quality and adapter trimming

```
mkdir trimmedData/

trim_galore --path_to_cutadapt /path/to/cutadapt \
    -o trimmedData/ \
    --paired sewage_sample_R1_001.fastq.gz sewage_sample_R2_001.fastq.gz > sewage_sample_trimming.output 2>&1
```

### Mapping to reference genome

```
minimap2 -t 8 \
    -ax sr NC_045512.fasta \
    trimmedData/sewage_sample_R1_001_val_1.fq.gz \
    trimmedData/sewage_sample_R2_001_val_2.fq.gz > sewage_sample.sam
```

### BAM file management

```
samtools view -@ 8 -Sb sewage_sample.sam > sewage_sample.bam
samtools sort -@ 8 -n sewage_sample.bam -o sewage_sample_sorted.bam

samtools view -@ 8 -u -f 1 -F 12 sewage_sample_sorted.bam > sewage_sample_mapped.bam
samtools sort -@ 8 sewage_sample_mapped.bam -o sewage_sample_sorted_mapped.bam
```

### Primer trimming

```
samtools index sewage_sample_sorted_mapped.bam

ivar trim \
    -b nCoV-2019.primer.bed \
    -q 15 \
    -m 200 \
    -s 4 \
    -i sewage_sample_sorted_mapped.bam \
    -p sewage_sample_sorted_mapped_trimmed.bam > sewage_sample_ivar.out
```

### Re-mapping to reference genome

```
samtools sort -@ 8 -n sewage_sample_sorted_mapped_trimmed.bam -o sewage_sample_sorted_mapped_trimmed_f.bam
bedtools bamtofastq -i sewage_sample_sorted_mapped_trimmed_f.bam -fq sewage_sample_R1_mapped.fastq -fq2 sewage_sample_R2_mapped.fastq

minimap2 -t 8 -ax sr NC_045512.fasta sewage_sample_R1_mapped.fastq sewage_sample_R2_mapped.fastq > sewage_sample_f.sam

samtools view -@ 8 -Sb sewage_sample_f.sam > sewage_sample_f.bam 
samtools sort -@ 8 -o sewage_sample_sorted_f.bam sewage_sample_f.bam 
```

### Filtering out duplicates

```
java -jar picard.jar MarkDuplicates I=sewage_sample_sorted_f.bam O=sewage_sample_sorted_uniq.bam M=sewage_sample_dup_metrics.txt
```

### Produce and annotate VCF file

```
freebayes \
    -f NC_045512.fasta \
    -F 0.01 \
    -C 1 \
    --pooled-continuous sewage_sample_sorted_uniq.bam > sewage_sample_freebayes.vcf

snpEff ann NC_045512.2 sewage_sample_freebayes.vcf > sewage_sample_freebayes_ann.vcf
```

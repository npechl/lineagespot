# Detecting SARS-CoV-2 lineages and mutational load in municipal wastewater; raw data analysis


## Steps


### Quality and adapter trimming

```
mkdir trimmedData/

trim_galore --path_to_cutadapt /path/to/cutadapt \
    -o trimmedData/ \
    --paired SampleA_R1_001.fastq.gz SampleA_R2_001.fastq.gz > SampleA_trimming.output 2>&1
```

### Mapping to reference genome

```
minimap2 -t 8 \
    -ax sr NC_045512.fasta \
    trimmedData/SampleA_R1_001_val_1.fq.gz \
    trimmedData/SampleA_R2_001_val_2.fq.gz > SampleA.sam
```

### BAM file management

```
samtools view -@ 8 -Sb SampleA.sam > SampleA.bam
samtools sort -@ 8 -n SampleA.bam -o SampleA_sorted.bam

samtools view -@ 8 -u -f 1 -F 12 SampleA_sorted.bam > SampleA_mapped.bam
samtools sort -@ 8 SampleA_mapped.bam -o SampleA_sorted_mapped.bam
```

### Primer trimming

```
samtools index SampleA_sorted_mapped.bam

ivar trim \
    -b nCoV-2019.primer.bed \
    -q 15 \
    -m 200 \
    -s 4 \
    -i SampleA_sorted_mapped.bam \
    -p SampleA_sorted_mapped_trimmed.bam > SampleA_ivar.out
```

### Re-mapping to reference genome

```
samtools sort -@ 8 -n SampleA_sorted_mapped_trimmed.bam -o SampleA_sorted_mapped_trimmed_f.bam
bedtools bamtofastq -i SampleA_sorted_mapped_trimmed_f.bam -fq SampleA_R1_mapped.fastq -fq2 SampleA_R2_mapped.fastq

minimap2 -t 8 -ax sr NC_045512.fasta SampleA_R1_mapped.fastq SampleA_R2_mapped.fastq > SampleA_f.sam

samtools view -@ 8 -Sb SampleA_f.sam > SampleA_f.bam 
samtools sort -@ 8 -o SampleA_sorted_f.bam SampleA_f.bam 
```

### Filtering out duplicates

```
java -jar picard.jar MarkDuplicates I=SampleA_sorted_f.bam O=SampleA_sorted_uniq.bam M=SampleA_dup_metrics.txt
```

### Produce and annotate VCF file

```
freebayes \
    -f NC_045512.fasta \
    -F 0.01 \
    -C 1 \
    --pooled-continuous SampleA_sorted_uniq.bam > SampleA_freebayes.vcf

snpEff ann NC_045512.2 SampleA_freebayes.vcf > SampleA_freebayes_ann.vcf
```
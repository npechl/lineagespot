# Commands used to generate VCF files

# mpileup tool

bcftools mpileup \
--fasta-ref /work/nikospech/SARS-CoV-2/ref/NC_045512.fasta \
Sewage_CoV19_L3_S1_L001_sorted_uniq.bam \
-o Sewage_CoV19_L3_S1_L001_mpileup.vcf

# freebayes tool

freebayes \
-f /work/nikospech/SARS-CoV-2/ref/NC_045512.fasta \
-F 0.5 \
-C 1 \
--pooled-continuous CoV19_L1_S1_sorted_uniq.bam > CoV19_L1_S1_freebayes_50.vcf



java -jar /work/precmed/varCall/picard.jar AddOrReplaceReadGroups \
I=Sewage_CoV19_L3_S1_L001_sorted_uniq.bam \
O=Sewage_CoV19_L3_S1_L001_sorted_uniq_rg.bam \
RGID=1 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=Sewage_CoV19_L3_S1_L001
	
samtools index Sewage_CoV19_L3_S1_L001_sorted_uniq_rg.bam 

/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java \
-jar /mnt/home/bio_tmp/apps/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar Mutect2 \
--reference /work/nikospech/SARS-CoV-2/ref/NC_045512.fasta \
--input Sewage_CoV19_L3_S1_L001_sorted_uniq_rg.bam \
--output Sewage_CoV19_L3_S1_L001_gatk_mutect2.vcf 

source activate covid19

snpEff ann NC_045512.2 CoV19_L1_S1_freebayes_10.vcf > CoV19_L1_S1_freebayes_10_ann.vcf

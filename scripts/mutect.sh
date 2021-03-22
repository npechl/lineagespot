/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java \
-jar /work/precmed/varCall/picard.jar AddOrReplaceReadGroups \
I=/work/nikospech/SARS-CoV-2/run-2021-02-01/CoV19_L1_S1_sorted_uniq.bam \
O=CoV19_L1_S1_sorted_uniq_rg.bam \
RGID=1 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=CoV19_L1_S1

samtools index CoV19_L1_S1_sorted_uniq_rg.bam 

/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java \
-jar /mnt/home/bio_tmp/apps/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar Mutect2 \
--reference /work/nikospech/SARS-CoV-2/ref/NC_045512.fasta \
--input CoV19_L1_S1_sorted_uniq_rg.bam \
--tumor-sample CoV19_L1_S1 \
--output CoV19_L1_S1_GATK_4_1_0_0_variants.vcf > CoV19_L1_S1.Mutect2_4_1_0_0.out 2>&1 
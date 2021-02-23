# Commands used to generate VCF files

# mpileup tool

bcftools mpileup \
--fasta-ref /work/nikospech/SARS-CoV-2/ref/NC_045512.fasta \
L1_S22_L001_sorted_uniq.bam \
-o L1_S22_L001_mpileup.vcf

# freebayes tool

freebayes \
-f /work/nikospech/SARS-CoV-2/ref/NC_045512.fasta \
-F 0.01 \
-C 1 \
--pooled-continuous Sewage-L2_S10_L001_sorted_uniq.bam > Sewage-L2_S10_L001_freebayes.vcf


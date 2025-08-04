# To create the modified GTF annotation file:
# 1 - Download Mus_musculus.GRCm38.102.chr.gtf.gz from Ensembl (GRCm38 and not GRCm39)
# 2 - Decompress (unzip) GTF file using the code below
gunzip Mus_musculus.GRCm38.102.chr.gtf.gz
# 3 - Remove lines corresponding to Gm45623 using the code below
grep 'gene_name "Gm45623"' Mus_musculus.GRCm38.102.chr.gtf -v > Mus_musculus.GRCm38.102.chr_irlab_annotation_ensembl102-20240302.smim32refseq_gtflines.gtf
# 4 - Add annotation corresponding to Smim32 using the code below
cat irlab_annotation_ensembl102-20240302.smim32refseq_gtflines.gtf >> Mus_musculus.GRCm38.102.chr_irlab_annotation_ensembl102-20240302.smim32refseq_gtflines.gtf

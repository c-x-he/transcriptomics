# Cellranger code is run in terminal using server (node 10)

# STEP 1: Download cell ranger and reference
curl -o cellranger-7.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.1.tar.gz?Expires=1669713504&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Njk3MTM1MDR9fX1dfQ__&Signature=eMdkKvu4Xq5GDMB87XpKkfeyvcxmXaIgPgoLuJM93hhPjkk0ejYW0aS1CLHUG-NUbKzDazFjC-7i5w36DdbDjsVa5yheDiJiJPFzDLBkv7RqfkYTyp00oU8a~iAF-2QkDfmy2ceUgEQbQCllgEHNYdmnpXmjScyW1UM2gzjYoOk~p2fBoWljqeg-qSGpAetftPJ7sxx7h4aSikpqKklkieZM6K-88rkakl48uEs8yaudG3fAHPJSHu8BW2ak7sTMXgiHL9RDnBbZjRooz95Sn17ZT1GQ~svmOGPZRuCfBQikWhkL4M-ka1Wm0y~Ij4KQ20iJP-t5W~LHAee4Z7UsPw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

## human reference
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

## unzip gz file
cd /opt

## unpack cell ranger
tar -xzvf cellranger-7.0.1.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz


# STEP 2: Set up cellranger
## Export PATH
export PATH=/opt/cellranger-7.0.1:$PATH

## Run cell ranger testrun
cellranger testrun --id=check_install

##recall cellrager
pwd cellranger-7.0.1
/home/c_x_he/cellranger-7.0.1

## Run cellranger mkref
cellranger mkref --genome=GRCh38-3.0.0_premrna --fasta=refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa --genes=GRCh38-3.0.0.premrna.gtf

cellranger --transcriptome=/home/c_x_he/GRCh38-3.0.0_premrna


# STEP 3: Use Cellranger
#### usage
#cellranger count 
#--id= name of the directory that will be created with the related archives
#--fastqs=path to the fastq files
#--sample= should have the exact same name as the prefix in your fastq.
#--transcriptome=path to the transcriptome

## Example code for how to set up  command
## Run a line of cellranger count for each sample you want to process
## Each sample should take 2-3 hours to process. To make things easier, you can run everything in a screen by typing screen -S screen_name and inputting everything into the screen so it keeps running even when your computer sleeps. To exit a screen, hit control-A-D, and to return to the screen, type screen -r screen_name

cellranger count --id=run_countpremRNA_C226_C01 --fastqs=/home/Camila/Mapping/C226_C01 --sample=226_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4






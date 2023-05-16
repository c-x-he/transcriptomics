
##### Download cell ranger #####
curl -o cellranger-7.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.1.tar.gz?Expires=1669713504&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Njk3MTM1MDR9fX1dfQ__&Signature=eMdkKvu4Xq5GDMB87XpKkfeyvcxmXaIgPgoLuJM93hhPjkk0ejYW0aS1CLHUG-NUbKzDazFjC-7i5w36DdbDjsVa5yheDiJiJPFzDLBkv7RqfkYTyp00oU8a~iAF-2QkDfmy2ceUgEQbQCllgEHNYdmnpXmjScyW1UM2gzjYoOk~p2fBoWljqeg-qSGpAetftPJ7sxx7h4aSikpqKklkieZM6K-88rkakl48uEs8yaudG3fAHPJSHu8BW2ak7sTMXgiHL9RDnBbZjRooz95Sn17ZT1GQ~svmOGPZRuCfBQikWhkL4M-ka1Wm0y~Ij4KQ20iJP-t5W~LHAee4Z7UsPw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

#human reference
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz


#unzip gz file
cd /opt

##### unpack cell ranger #####
tar -xzvf cellranger-7.0.1.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

##### Export PATH #####
export PATH=/opt/cellranger-7.0.1:$PATH

##### Run cell ranger testrun #####
cellranger testrun --id=check_install


#recall cellrager
pwd cellranger-7.0.1
/home/c_x_he/cellranger-7.0.1

##### Create a premRNA GTF ##### ### You don't need to run this step, since the premRNA is already created. 
# (1)
#awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf > GRCh38-3.0.0.premrna.gtf

# (2) Run cellranger mkref
#cellranger mkref --genome=GRCh38-3.0.0_premrna --fasta=refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa --genes=/home/Camila/refdatacellranger-GRCh38-3.0.0/genes/GRCh38-3.0.0.premrna.gtf
#cellranger mkref --genome=GRCh38-3.0.0_premrna --fasta=refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa --genes=GRCh38-3.0.0.premrna.gtf

cd /home/c_x_he/SingleSplice/

cellranger mkref --genome=GRCh38-3.0.0_premrna_MAPT4 \
  --fasta=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna/fasta/genome.fa \
  --genes=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna/genes/genes.gtf 

cellranger --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

##### Download results from sequencing ##### 

#replace the number project after #Data/ with the project you want to download

wget -r -nH -nc -R index.html* http://slimsdata.genomecenter.ucdavis.edu/Data/tkiw9nfdnj


##### Cellranger count #####

### cellranger will be at your working directory

#### usage
cellranger count 
--id= name of the directory that will be created with the related archives
--fastqs=path to the fastq files
--sample= should have the exact same name as the prefix in your fastq.
--transcriptome=path to the transcriptome

####

cellranger count --id=run_countpremRNA_C226_C01 --fastqs=/home/Camila/Mapping/C226_C01 --sample=226_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C240_C01 --fastqs=/home/Camila/Mapping/C240_C01 --sample=240_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C262_C01 --fastqs=/home/Camila/Mapping/C262_C01 --sample=262_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C291_C01 --fastqs=/home/Camila/Mapping/C291_C01 --sample=291_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C327_C01 --fastqs=/home/Camila/Mapping/C327_C01_Nova --sample=C327_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C342_C01 --fastqs=/home/Camila/Mapping/C342_C01 --sample=342_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C350_C01 --fastqs=/home/Camila/Mapping/C350_C01_All --sample=350_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C364_C01 --fastqs=/home/Camila/Mapping/C364_C01 --sample=364_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C114_C01 --fastqs=/home/Camila/Mapping/C114_C01 --sample=114_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C135_C01 --fastqs=/home/Camila/Mapping/C135_C01 --sample=135_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C140_C01 --fastqs=/home/Camila/Mapping/C140_C01 --sample=140_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C142_C01 --fastqs=/home/Camila/Mapping/C142_C01 --sample=142_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C361_C01 --fastqs=/home/Camila/Mapping/C361_C01 --sample=c361_c1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C363_C01 --fastqs=/home/Camila/Mapping/C363_C01 --sample=363_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C379_C01 --fastqs=/home/Camila/Mapping/C379_C01 --sample=379_C1 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4


cellranger count --id=run_countpremRNA_C381_C01 --fastqs=/home/Camila/Mapping/C381_C01 --sample=C381 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C323_C01 --fastqs=/home/Camila/Mapping/C323_C01 --sample=C323 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4


cellranger count --id=run_countpremRNA_C255_C01 --fastqs=/home/Camila/Mapping/C255_C01 --sample=C255 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4


cellranger count --id=run_countpremRNA_C210_C01 --fastqs=/home/Camila/Mapping/C210_C01 --sample=C210 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4



cellranger count --id=run_countpremRNA_C269_C01 --fastqs=/home/Camila/Mapping/C269_C01 --sample=C269 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C99_C01 --fastqs=/home/Camila/Mapping/C99_C01 --sample=C99 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4



cellranger count --id=run_countpremRNA_C181_C01 --fastqs=/home/Camila/Mapping/C181_C01 --sample=C181 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4

cellranger count --id=run_countpremRNA_C277_C01 --fastqs=/home/Camila/Mapping/C277_C01 --sample=C277 --transcriptome=/home/c_x_he/SingleSplice/GRCh38-3.0.0_premrna_MAPT4



























### for each sample you will run a code from the above. Each sample takes about 2-3 hours to run, it depends on the memory available on the server, so what do I do, I write the script for all the samples and paste it to run at once on the server, so when first sample is finished, it automatically switches to the next sample.

####Good luck!!!!###






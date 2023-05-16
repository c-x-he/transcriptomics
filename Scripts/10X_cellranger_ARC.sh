cp -r /home/c_x_he/cellranger-arc-2.0.2 /home/(your folder)/cellranger-arc-2.0.2

cp -r /home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 /home/(your folder)/refdata-cellranger-arc-GRCh38-2020-A-2.0.0


#download cellranger_arc and human reference
#curl -o cellranger-arc-2.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1669791323&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hcmMvY2VsbHJhbmdlci1hcmMtMi4wLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjY5NzkxMzIzfX19XX0_&Signature=iPvBM2TPuzVEusaBt4oICdzvxUimGqRqvxaCH5GTP1F~K16utIHfntTo160vBfw962FUU52dYOGtj7dUSWrRHqCpD4tSjAuOFKiu5ylYNkry4AKN9dLFvAh5TyjYTxVpo08Tm4HyhbYOqGfufyJeZUuVeE526aNUs3k9k18~WPNydz2qBAL37lEIuiDsuC~34aIenAYgxpPC~y98EP-MwXzqcWnHY~4AgP~UoaihL0Y4C8VrSzDtqlM2G0RYRM-ChaQAhb5aAseywHpjY3QzevjKvNO6-USpLoLpttqKDW2rX3wdtP-fbQSWbZHbXlx-XCFpykJCNz7TnK5175XLZg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

#tar -xzvf cellranger-arc-2.0.2.tar.gz
 
#curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

#tar -xzvf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

export PATH=/home/c_x_he/cellranger-arc-2.0.2:$PATH
 
#cellranger-arc sitecheck > sitecheck.txt
#cellranger-arc upload c_x_he@ucsb.edu sitecheck.txt

#cellranger-arc testrun --id=tiny

#make csv file using excel or something with 3 columns: fastqs, sample, and library_type
#fill it out like in: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/count
cd /home/c_x_he/multiome/outs

cellranger-arc count --id=364_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/364_library.csv \
                       --min-atac-count=1000 \
                       --min-gex-count=200

cellranger-arc count --id=321_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/321_library.csv \
                       --min-atac-count=700 \
                       --min-gex-count=2000


cellranger-arc count --id=299_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/299_library.csv \
                       --min-atac-count=800 \
                       --min-gex-count=300


cellranger-arc count --id=201_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/201_library.csv \
                       --min-atac-count=500 \
                       --min-gex-count=500


cellranger-arc count --id=160_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/160_library.csv \
                       --min-atac-count=500 \
                       --min-gex-count=500
                       

cellranger-arc count --id=147_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/147_library.csv \
                       --min-atac-count=700 \
                       --min-gex-count=2000


cellranger-arc count --id=118_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/118_library.csv \
                       --min-atac-count=500 \
                       --min-gex-count=1000



cellranger-arc count --id=109_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/109_library.csv \
                       --min-atac-count=200 \
                       --min-gex-count=500









###################################################
cellranger-arc count --id=364_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/364_library.csv \
                       --min-atac-count=2000 \
                       --min-gex-count=500
                       
                       
                       
cellranger-arc count --id=118_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/118_library.csv \
                       --min-atac-count=1000 \
                       --min-gex-count=1000


cellranger-arc count --id=160_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/160_library.csv \
                       --min-atac-count=1000 \
                       --min-gex-count=700
                       




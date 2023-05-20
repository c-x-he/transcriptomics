# STEP 1: SET UP ENVIORNMENT

## 1. download cellranger_arc and human reference
curl -o cellranger-arc-2.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1669791323&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hcmMvY2VsbHJhbmdlci1hcmMtMi4wLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjY5NzkxMzIzfX19XX0_&Signature=iPvBM2TPuzVEusaBt4oICdzvxUimGqRqvxaCH5GTP1F~K16utIHfntTo160vBfw962FUU52dYOGtj7dUSWrRHqCpD4tSjAuOFKiu5ylYNkry4AKN9dLFvAh5TyjYTxVpo08Tm4HyhbYOqGfufyJeZUuVeE526aNUs3k9k18~WPNydz2qBAL37lEIuiDsuC~34aIenAYgxpPC~y98EP-MwXzqcWnHY~4AgP~UoaihL0Y4C8VrSzDtqlM2G0RYRM-ChaQAhb5aAseywHpjY3QzevjKvNO6-USpLoLpttqKDW2rX3wdtP-fbQSWbZHbXlx-XCFpykJCNz7TnK5175XLZg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -xzvf cellranger-arc-2.0.2.tar.gz
 
curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

tar -xzvf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

export PATH=/home/c_x_he/cellranger-arc-2.0.2:$PATH

## 2. make csv file using excel or something with 3 columns: fastqs, sample, and library_type
## fill it out like in: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/count

# STEP 2: RUN CELLRANGER-ARC 
## 1. Specify working directory
cd /home/c_x_he/multiome/outs

## 2. Run cellranger-arc
## run this code for each sample you want processed, the basic structure of this command is 
cellranger-arc count --id=[sample ID]\
                       --reference=[path to GRCh38-2020-A-2.0.0] \
                       --libraries=[path to the csv file created in step 1.2] \
                       --min-atac-count= [optional, manually specifcy minimum number of atac counts to be considered a cell] \
                       --min-gex-count= [optional, manually specifcy minimum number of GEX counts to be considered a cell]
## here is an example
cellranger-arc count --id=364_C04\
                       --reference=/home/c_x_he/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/home/c_x_he/multiome/364_library.csv \
                       --min-atac-count=1000 \
                       --min-gex-count=200




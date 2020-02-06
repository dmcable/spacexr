For cloning the repository:
> git clone git@github.com:dmcable/slideseq.git
To update:
> git pull origin master

To get workspace:
> git checkout -b workspace
> git pull origin workspace


List of packages to install: dplyr, pals
List of packages to update: 
Install a package:
> R
> install.packages(“package”)
> (ctr-C to exit R)
> install.packages("toupdate") #updates the package


In order to run the cell type decomposition method, the first step is to process the single cell reference. Create a folder e.g. 'SingleCellReference' containing the following three files:
1) meta_data.csv # a CSV file (with 3 columns, with headers "barcode", "cluster", and "nUMI") containing the numeric cluster assignment for each cell
2) cell_type_dict.csv # a CSV file (with 2 columns, with headers "Cluster" and "Name") containing the mapping between numeric cluster ID and cluster name. If you want a cluster to be filtered out of the single cell reference, you can leave the cluster name blank. The cell types must not contain the character '/' or '-'
3) dge.csv # a DGE (barcodes by gene counts) CSV file in the standard 10x format

Run the script dgeToSeurat.R (changing appropriately mydir and refdir). This creates a Seurat object stored in SCRef.RDS in your data directory.
> Rscript dgeToSeurat.R

For the python script, you need to install tensor flow:
> pip install tensor flow

Change the PEERCTD.sh script to executable:
> chmod u+x PEERCTD.sh

Next, run the script to find the marker genes in the reference (changing the directory). If you have a small number of genes (e.g. <2500), then it is not necessary, and you could try using all the genes. Note: findMarkerGenes.R requires > 16G of memory (works on 64G).
> Rscript findMarkerGenes.R
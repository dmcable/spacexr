To run this on Broad server, use:
> ish -l os=RedHat7 -l h_vmem=32G
> use R-3.5
> use Python-3.6
Note: these versions don't need to be exact
You must install pandas, tensorflow, and yaml package in Python before using the RCTD package:
> pip install --upgrade --user tensorflow
> pip install --upgrade --user pandas
> pip install --user --upgrade pyyaml

For cloning the repository:
> git clone git@github.com:dmcable/RCTD.git

Make sure you are on branch dev:
> git status
> git checkout dev
> git status

Next, make sure that you have all of the packages Installed in R that are listed under Imports in RCTD/DESCRIPTION
Next, build the package:
> R CMD build RCTD
This will generate a file such as RCTD_0.1.0.tar.gz. Run R (e.g. 3.5) and install the package:
> install.packages("RCTD_0.1.0.tar.gz", repos = NULL)
Now, you should be able to use:
> library(RCTD)

To update:
> git pull origin dev

To get a sample Data directory:
> cp -r /broad/thechenlab/Dylan/slideseq/Cell Demixing/SampleData/Data YOUR_RCTD/

For an example on how to format for dgeToSeurat, see SampleData/Data/Reference/KidneyHumanReference
The Data folder should have 'Reference' and 'Slideseq' subdirectories. 

In order to run the cell type decomposition method, the first step is to process the single cell reference. Create a folder in 'Data/Reference' e.g. 'Data/Reference/SingleCellReference' containing the following three files:
1) meta_data.csv # a CSV file (with 3 columns, with headers "barcode", "cluster", and "nUMI") containing the numeric cluster assignment for each cell
2) cell_type_dict.csv # a CSV file (with 2 columns, with headers "Cluster" and "Name") containing the mapping between numeric cluster ID and cluster name. If you want a cluster to be filtered out of the single cell reference, you can leave the cluster name blank. The cell types must not contain the character '/' or '-'
3) dge.csv # a DGE (barcodes by gene counts) CSV file in the standard 10x format

Run the function dgeToSeurat(refdir) This creates a Seurat object stored in SCRef.RDS in your reference data directory:
> dgeToSeurat('Data/Reference/KidneyHumanReference')

Next, put the slideseq puck folder in your Data/Slideseq directory. The puck folder needs to contain the BeadLocationForR.csv file as well as MappedDgeForR.csv (or other named DGE, you can change the name with the puckfile option in config.yml).

VERY IMPORTANT: before running PEERCTD, make sure to edit the config.yml file. There are some initial parameters recommended for testing the package will run. Ultimately, you should change to the commented out parameters for the true run.

Without changing config.yml, you can run PEERCTD on the sample data:
> ./PEERCTD
After running, you should see the following output files in Data/Slideseq/YOURPUCK/results:
-all_cell_types.pdf, an image of all the cell types
-cell_type_calls.pdf, the individual cell types
-puck_labeled.RDS an RDS containing the puck labeled with cell types
-test_results.RDS an RDS containing the results of the cell type calling procedure.


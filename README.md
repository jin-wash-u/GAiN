# GAiN
GAiN generates large synthetic cohorts using small training sets of gene expression data from samples of two different phenotypes. It then performs differential expression (DE) analysis on the synthetic cohorts, returning a list of candidate DE genes.
   
## Prerequisites
Please make sure you have installed the following tools:
- [Python3](https://www.python.org/)
  - [pandas](https://pandas.pydata.org/)
  - [tensorflow](https://www.tensorflow.org/)
  - [keras](https://keras.io/)
  - [sklearn](https://scikit-learn.org/)
- [R](https://www.r-project.org/)
  - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

## Installation
Clone this repository to the desired location:

```
git clone https://github.com/jin-wash-u/GAiN.git
```

## Test installation
To test the installation, execute the following command in your terminal to display usage information about GAiN. 
```
$PATH_TO_GAIN/GAiN -h
```

## Parameters
Running GAiN with the -h option (or --help) will print a desciption of its optional and required input arguments. A description of each follows.
```
GAiN [-h] [--version] [-b C0 C1] [-e EPOCHS] [-a ALPHAGAN] [--minGE GE]
     [--minLFC LFC] [--numbOfNetworks NON] [--numNetworkCutoff NNC]
     [--deseq] [--save] [--synth] [--seed SEED] [-q] [-o OUTNAME]
     input.csv population.csv
```
### Optional arguments
- -h, --help  
     Prints the help menu, which contains an example of the function usage and abbreviated explanations of each of the options.
- --version  
     Prints version information
- -b C0 C1, --batchsizes C0 C1  
     Size of synthetic cohort to generate for each condition (default: 500 500)
- -e EPOCHS, --epochs EPOCHS  
     Number of epochs for training model (default: 500)
- -a ALPHA, --alphaGAN ALPHA  
     Significance threshold for reporting genes as DE between the synthetic groups (default: 0.05)
- --minGE GE  
     Minimum absolute difference in average expression between the synthetic groups for a gene to be reported as DE (default: 10)
- --minLFC LFC  
     Minimum log2 fold change in expression between the synthetic groups for a gene to be reported as DE (default: 1)
- --numbOfNetworks NON  
     Number of networks to use for bagging (default: 20)
- --numNetworkCutoff NNC  
     Number of networks gene must be significantly DE in (default: 20)
- --deseq  
     Use DESeq2 method for DE significance calculations (default: use edgeR)
- --save  
     Save trained models for later use
- --synth  
     Save synthetic expression tables
- --seed SEED    
     Optional seed for random sampling of the training sets
- -q, --quiet  
     Run in quiet mode, limiting program output
- -o OUTNAME, --outname OUTNAME  
     Prefix for output filenames (default: ./GAiN)

### Input files
- input.csv  
     Path to the training cohort gene expression table in CSV format.  A pair of large synthetic cohorts will be generated based on the samples in this table.  The first row must be any label string followed by the sample ids of the training set.  The second row must start with a label string (e.g. "Cancer ID") followed by 0 or 1 based on the phenotypic group of the sample for that column.  All subsequent rows must start with a gene label, followed by the expression level of that gene in the sample for that column.
- population.csv  
     Path to the population cohort gene expression table in CSV format.  Scale will be restored to the synthetic gene expression tables using these samples.  Follows the same format as input.csv, but without any phenotype group/Cancer ID row.  Note that only genes with entries in both tables will be modeled.

### Output
GAiN generates a CSV table, by default GAiN_DE_genes.csv, containing the list of differentially expressed genes between the two synthetic groups that passed all filters, together with the sum of each gene's rank across the NON networks and the number of networks in which it was significantly DE.
 
### Example
An example input expression CSV file is included with GAiN to demonstrate how to run the tool. It can be run as follows:

```
cd $PATH_TO_DANSR
GAiN \
        -o test \
	--seed 42 \
	example/example_input.csv \
	example/example_population.csv
```
The results of this execution can be compared with the results file example/example_result_DE_genes.csv.

The gene expression values in example_input.csv consist of TMM-normalized counts from 10 primary tumors of luminal B subtype (condition 0) and 10 primary tumors of triple negative subtype (condition 1) sourced from TCGA-BRCA.  The population cohort is TMM normalized expression of all primary tumors in TCGA-BRCA.

## GAiN Docker Instructions
A docker image for GAiN has been created and tested on Linux and Mac. To run GAiN using this method, you need to have [Docker](https://docs.docker.com/) installed on your machine. 

### Docker Installation
* Ubuntu: follow [the instructions](https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/) to get Docker CE for Ubuntu.
* Mac: follow [the instructions](https://store.docker.com/editions/community/docker-ce-desktop-mac) to install [the stable version of Docker CE](https://download.docker.com/mac/stable/Docker.dmg) on Mac.
<!--- 
* Windows: follow [the instructions](https://docs.docker.com/toolbox/toolbox_install_windows/) to install [Docker Toolbox](https://download.docker.com/win/stable/DockerToolbox.exe) on Windows. 
-->
 
To obtain the latest docker image, run the following on your command line:
 
```
docker pull mjinkm/gain
```
To test the image, run the following command which shows the usage of this tool:
```
docker run mjinkm/gain GAiN -h
```
### Example 
To run GAiN using docker on the example provided, use the following command:
```
docker run -v $PATH_TO_OUTPUT:/gain_out mjinkm/gain GAiN \
        -o /gain_out/test \
	--seed 42 \
	/opt/GAiN/example/example_input.csv \
	/opt/GAiN/example/example_population.csv
```
where `PATH_TO_OUTPUT` is the desired local output directory.

# GAiN
GAiN generates large synthetic cohorts using sparse training sets of gene expression data from samples of two different phenotypes. It then performs differential expression (DE) analysis on the synthetic cohorts, returning a list of candidate DE genes.
   
## Prerequisites
Please make sure you have installed the following tools:
- [Python3](https://www.python.org/)
  - [pandas](https://pandas.pydata.org/)
  - [tensorflow](https://www.tensorflow.org/)
  - [keras](https://keras.io/)

## Installation
Clone this repository to the desired location:

```
git clone https://github.com/mjinkm/GAiN.git
```

## Test installation
To test the installation, execute the following command in your terminal to display usage information about GAiN. 
```
$PATH_TO_GAIN/GAiN -h
```

## Parameters
Running HPV-EM with the -h option (or --help) will print a desciption of its optional and required input arguments. A description of each follows.
```
GAiN [-h] [-s C0 C1] [-b C0 C1] [-e EPOCHS] [-a ALPHA] [--minGE GE] [--minLFC LFC]
     [--seed SEED] [-q] [-o OUTNAME] input.csv
```
### Optional arguments
- -h, --help
     Prints the help menu, which contains an example of the function usage and abbreviated explanations of each of the options.
- -s C0 C1, --samplenums C0 C1
     Number of samples from each condition to use in training GAN (set to -1 to use all from a condition (default: -1 -1)
- -b C0 C1, --batchsizes C0 C1
     Size of synthetic cohort to generate for each condition (default: 1000 1000)
- -e EPOCHS, --epochs EPOCHS
     Number of epochs for training model (default: 80)
- -a ALPHA, --alphaGAN ALPHA
     Significance threshold for reporting genes as DE between the synthetic groups (default: 1.25e-8/(total number of genes tested)
- --minGE GE  
     Minimum absolute difference in average expression between the synthetic groups for a gene to be reported as DE (default: 10)
- --minLFC LFC
     Minimum log2 fold change in expression between the synthetic groups for a gene to be reported as DE (default: 1)
- --seed SEED
     Optional seed for random sampling of the training sets
- -q, --quiet
     Run in quiet mode, limiting program output
- -o OUTNAME, --outname OUTNAME
     Prefix for output filenames (default: ./GAiN)

### Input file
- input.csv
     Path to a CSV table.  The first row must be any label string followed by the sample ids of the training set.  The second row must start with "Cancer ID" followed by 1 or 0 based on qhich phenotypic group of the sample for that column.  All subsequent rows must start with a gene label, followed by the expression level of that gene in the sample for that column.

### Output
GAiN generates a CSV table, by default GAiN_DE_genes.csv, containing the list of differentially expressed genes between the two synthetic groups that passed the significance, absolute difference and LFC filters, along with the log2 fold change and nominal p-value for each.
 
### Example
An example input expression CSV file is included with GAiN to demonstrate how to run the tool. It can be run as follows:

```
cd $PATH_TO_DANSR
GAiN \
        -s 10 10 \
        -o test \
	example/example_input.csv
```
The results of this execution can be compared with the results file example/example_result_DE_genes.csv.

The gene expression values in example_input.csv sourced from GSE22260 consisting of 19 prostate tumor (condition 0) and 10 adjacent normal (condition 1) samples.

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
        -s 10 10 \
        -o /gain_out/test \
	/opt/GAiN/example/example_input.csv
```
where `PATH_TO_OUTPUT` is the desired local output directory.

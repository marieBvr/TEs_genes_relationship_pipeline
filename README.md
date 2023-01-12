# TEGRiP: Transposable Elements Genes RelationshIps Pipeline

TEGRiP is recommended by PCI Genomics [10.24072/pci.genomics.100010](https://doi.org/10.24072/pci.genomics.100010) so please cite [10.1101/2021.02.25.432867](https://doi.org/10.1101/2021.02.25.432867) if you use it.

## Context

The main goal of this project is to find the positional relationships between Transposable Elements (TEs) and genes along genome. 

The Apricot genome annotation has been used to validate our strategy (raw data available on ENA [PRJEB42606](https://www.ebi.ac.uk/ena/browser/view/PRJEB42606)). This pipeline can be used with custom TE annotation as well as de novo assembled genome of any kind of species.

You can learn more about this subject in the Contribution Section.

## Introduction 

Transposable elements are DNA fragment capable of moving from one place to another troughout the genome  via a mecanism called transposition. 
There are different category/class of transposon. In this project we are going to focus on LTR (long terminal repeat). Learn more by [clicking here ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874221/)

## Requirement

To use this programme please clone the repository. Make sure that Python3 and Rstudio are installed on your machine. 
If it's not the case, check out the following links :
- python 3 : [link](https://www.python.org/downloads/)
- Rstudio/R >= 4.0.2 : [link](https://rstudio.com/products/rstudio/download/)

# Example 

Before running the programme for your own data, please use the testing data to check that everything works.

In the `data` folder, there are 2 files, one with the data regarding each gene and one regarding each transposon displayed in the diagram above. 
The python scripts are going to find for each transposable element the nearest gene before and after it. The script will also look for overlapping. 
The result file of this test data is located in the `result` folder. 

# Usage
-----------------------
## Python Scripts 
### Input files
The programme takes in consideration the columns name so before using them please verify that each files have those columns : 

* Genome files

chromosome | source | feature | start | end | score | strand | phase | ID | Attributes | 
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- 

* TE files

Chromosome | Length_Chr | Type | X | Start | End | X | X | Strand | X | Attribute | X | Class | TE_name | X |
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- |--- |--- |--- 

* LTR files

species | ID | dfam_target_name | X | X | X | X | chromosome | start | end | strand | X | annotation | X | X | score |
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- | --- | --- | --- |---

### Pipeline

Scripts are time optimized using multiprocessing. 
Multiprocessing is a system that use multiple central processing units (CPUs) making the scripts run faster.

To run the script type the following line by replacing each file name by the real name.

To analyze general TE:

```bash
python3 Multiprocessing/Create_Data_TE_multipro_reformatted.py \
	-g data/Gene_testing_data.tsv \
	-te data/Transposon_testing_data.tsv \
	-o result/output_TE.tsv
```

To analyze LTR:

```bash
python3 Multiprocessing/Create_Data_LTR_multiprocessing_reformatted.py \
	-g data/Gene_testing_data.tsv \
	-te data/LTR_testing_data.tsv \
	-o result/output_LTR.tsv
```

The script will take each file and extract all the data and put them in lists of dictionaries. 
Then for each TE, it will check the nearest gene whether it's subset,superset, upstream/downstream or whether it's an upstream/downstream overlap.

### Output file

The script create a new output file (.tsv) which will be used to make the statistical analysis.

### Post analysis
In the first analysis, the strand (+ or -) of genes and TEs are considered.
It is possible to take it into accompt by running:
```bash
python3 Multiprocessing/Post_TE_multipro_reformatted.py -i result/output_TE.tsv -o result/output_post_TE.tsv
```

-----------------------
## R scripts

There are four R scripts allowing to report different kinds of information:
 - Count the number of TEs
 - Count the number of TEs with associated genes within a certain distance
 - Overlap statistics, which show how many TEs have an overlap with gene, both upstream and downstream.
 - Distance statistics, which show the number of each TE superfamily overlap with the closest gene (-i -1 -x 0) or within the distance of 0-500 bp (-i 0 -x 500), 500-1000 bp (-i 501 -x 1000), 1000-2000 bp (-i 1001 -x 2000) and more than 2000 bp. The subsets and supersets are not included in these counts. The input file is the result file obtained from the Python script.

The input file is the result file obtained with the Python script.

```bash
Rscript Rscript/number_te.r \
	-f result/output_TE.tsv \
	-o result/count_TE_transposons.pdf

Rscript Rscript/count_gene_associated_te.r \
	-f result/output_TE.tsv \
	-o result/count_TE_transposons.pdf \
	-i 0 \
	-x 2000

Rscript Rscript/Overlap_counting.r \
	-f result/output_TE.tsv \
	-p result/overlap_TE_results.pdf \
	-o result/overlap_TE_results.csv

Rscript Rscript/Distance_counting.r \
	-f result/output_TE.tsv \
	-p result/distance_TE_results.pdf \
	-o result/distance_TE_results.csv
```

Please note that the **graph's legend will also need to be change** according to the file and species.

# Contribution
This programme has been developped by [Caroline Meguerditchian](mailto:caroline.meguerditchian@etu.u-bordeaux.fr) and [Ayse Ergun](mailto:aergun@u-bordeaux.fr) 
under the supervision of [Marie Lefebvre](mailto:marie.lefebvre@inrae.fr) and [Quynh Trang-Bui](mailto:quynh-trang.bui@inrae.fr).





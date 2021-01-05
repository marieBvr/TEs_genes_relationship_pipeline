# TE_abricot

## Context

The main goal of this project is to find the relationship between Transposable Elements (TEs) and genes along genome. 
Then focus on the Long Terminal Repeat (LTR) that are type of TE. 

As an example we used the Apricot genome and compare TEs for 4 different species of Apricot.

You can learn more about this subject in the Contribution Section.

## Introduction 

Transposable elements are DNA fragment capable of moving from one place to another troughout the genome  via a mecanism called transposition. 
There are different category/class of transposon. In this project we are going to focus on LTR (long terminal repeat). Learn more by [clicking here ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874221/)

## Requirement

To use this programme please clone the repository. Make sure that Python3 and Rstudio are installed on your machine. 
If it's not the case, check out the following links :
- python 3 : [link](https://www.python.org/downloads/)
- Rstudio/R : [link](https://rstudio.com/products/rstudio/download/)

An access to a terminal is also required.

# Example 

Before running the programme for your own data, please use the testing data to check that everything works.

In the `Testing_data` folder, there are 2 files, one with the data regarding each gene and one regarding each transposon displayed in the diagram above. 
The python scripts are going to find for each transposable element the nearest gene before and after it. The script will also look for overlapping. 
The result file of this test data is located in the `Testing_result` folder. 

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

```
$ python3 Multiprocessing/Create_Data_multipro.py -g Testing_data/Gene_testing_data.tsv -te Testing_data/Transposon_testing_data.tsv -o Resting_result/ResultFile_TE.tsv

```

To analyze LTR:

```
$ python3 Multiprocessing/Create_Data_LTR_multiprocessing.py  -g [file_path]/gene.tsv  -te [file_path]/transposon.tsv -o [file_path]/output_LTR.tsv

```

The script will take each file and extract all the data and put them in lists of dictionaries. 
Then for each TE, it will check the nearest gene whether it's subset,superset, upstream/downstream or whether it's an upstream/downstream overlap.

### Output file

The script create a new output file (.tsv) which will be used to make the statistical analysis.

-----------------------
## R scripts
After opening the .r files with Rstudio, make sure to verify and modify the following line with the right file. 
The Input file is the result file obtained with the python script.

```
result_file = read.csv(file = '[file_path]/result.tsv', sep = '\t', header = TRUE)
```
Each .r script will give you a new output file (counter) as well as a graph. The title and small explication in the file will help to understand better.

Please note that the graph's legend will also need to be change according to the file and abricot species.

# Contribution
This programme has been developped by [Caroline Meguerditchian](caroline.meguerditchian@etu.u-bordeaux.fr) and [Ayse Ergun](aergun@u-bordeaux.fr) 
under the supervision of [Marie Lefebvre](marie.lefebvre@inrae.fr) and [Quynh Trang-Bui](quynh-trang.bui@inrae.fr).. 





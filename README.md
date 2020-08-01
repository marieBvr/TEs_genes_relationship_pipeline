# TE_abricot

## Context
This is the repository of our summer internship 2020. We are two students from University of Bordeaux. We had the chance of doing our internship in INRAE (french National Institute of Agricultural Research) alongside the Virology Team.

The main goal of this project is to find the relationship between Transposable Element (TE) and the apricot gene. Then focus on the Long Terminal Repeat (LTR) that are type of TE. And finally compare those element for 4 species of apricot.

You can learn more about this subject in the Contribution Section.

## Introduction 
Transposable element are DNA fragment capable of moving from one place to another troughout the genome  via a mecanism called transposition. 
There are different category/class of transposon. In this project we are going to focus on LTR (long terminal repeat). Learn more by [clicking here ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874221/)

## Requirement

To use these scripts please download the folder. Make sure that python3 and Rstudio are installed in your machine. If it's not the case, check out the following links :
- python 3 : [link](https://www.python.org/downloads/)
- Rstudio : [link](https://rstudio.com/products/rstudio/download/)

An access to a terminal is also required.

# Example 
![alt text](https://raw.githubusercontent.com/Ayse1006/TE_abricot/master/Testing_data/diagram_gene_te.jpg)

In the Testing_data folder, there are 2 files, one with the data regarding each gene and one regarding each transposon displayed in the diagram above. The python scripts are going to find for each transposable element the nearest gene before and after it. The script will also look for overlapping. The result file of this test data is located in the Testing_result folder. 

# Usage
-----------------------
## Python Scripts 
### Input files
The scripts takes in consideration the columns name so before using them please verify that each files haves thise columns : 

* Genome files

chromosome | source | feature | start | end | score | strand | phase | ID | Attributes | 
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- 

* TE files

Chromosome | Length_Chr | Type | match/match_part_part | Start | End | Length | Frame | Attribute | Code | Class | TE_name | TE_Status |
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- |--- 

* LTR files

species | ID | dfam_target_name | ltr_similarity | similarity | protein_domain | orfs | chromosome | start | end | strand | width | annotation | pred_tool | frame | score |
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- | --- | --- | --- |---

### Pipeline

To run the script just type the following line by replacing each file name by the real name.

Please use the Creat_Data_multipro.py if working with general TE and Create_Data_LTR_multiprocessing.py if working with LTR from the Multiprocessing folder. Multiprocessing is a system that use multiple central processing units (CPUs) making the scripts run faster.

```
$ python3 [script_path]/script.py  [file_path]/gene.tsv  [file_path]/transposon.tsv

```
The script will take each file and extract all the data and put them in lists of dictionaries. Then for each TE, it will check the nearest gene whether it's subset,superset, upstream/downstream or whether it's an upstream/downstream overlap.

### Output file

The script will creat a new output file(Testing_result folder) which will be used to make graphs.

-----------------------
## R scripts
After opening the .r files with Rstudio, make sure to verify and modify the following line with the right file. The Input file is the result file obtained with the python script.

```
result_file = read.csv(file = '[file_path]/result.tsv', sep = '\t', header = TRUE)
```
Each .r script will give you a new output file(counter) as well as a graph. The title and small explication in the file  will help to understand better.

# Contribution
You can find our reports on clicking on the following links.
Ayse's report :  https://fr.overleaf.com/read/qwygcwcjmxdd
Caroline's report :  https://www.overleaf.com/5677149684vhkkmwryfwjp






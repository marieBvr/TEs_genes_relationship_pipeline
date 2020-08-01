# TE_abricot

## Context
This is the repository of our summer internship 2020. We are two students from University of Bordeaux. We had the chance of doing our internship in INRAE (french National Institute of Agricultural Research) alongside the Virology Team.

The main goal of this project is to find the relationship between Transposable Element (TE) and the apricot gene. Then focus on the Long Terminal Repeat (LTR) that are type of TE. And finally compare those element for 4 species of apricot.

You can learn more about this subject in the Contribution Section.

## Introduction 
Transposable element are DNA fragment capable of moving from one place to another troughout the genome  via a mecanism called transposition. 
There are different category/class of transposon. In this project we are going to focus on LTR (long terminal repeat). You can learn more by [clicking here ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874221/)

## Requirement
In order to you these script you will have to download  our repository.
# Example 
![alt text](https://raw.githubusercontent.com/Ayse1006/TE_abricot/master/Testing_data/diagram_gene_te.jpg)

# Usage
```
$ python Noms_du_script gene_data te/ltr_data 
```
## Input files

The scripts takes in consideration the columns name so before using them you should verify that each files haves thise columns : 

* Genome files

chromosome | source | feature | start | end | score | strand | phase | ID | Attributes | 
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- 

* TE files

Chromosome | Length_Chr | Type | match/match_part_part | Start | End | Length | Frame | Attribute | Code | Class | TE_name | TE_Status |
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- |--- 

* LTR files

species | ID | dfam_target_name | ltr_similarity | similarity | protein_domain | orfs | chromosome | start | end | strand | width | annotation | pred_tool | frame | score |
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- | --- | --- | --- |---

## Output file
## R 
### Input files
### Output files

## Basic statistic ? (maybe ) 


# Contribution
You can find our reports on clicking on the following links.
Ayse's report :  https://fr.overleaf.com/read/qwygcwcjmxdd
Caroline's report :  https://www.overleaf.com/5677149684vhkkmwryfwjp






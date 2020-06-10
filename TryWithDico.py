#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Ayse Ergun
# Script For TE project

#==============================================================================
#                           importations
#==============================================================================
import sys
import csv
ListOfTE=[]
ListOfGene=[]

#==============================================================================
#                            functions
#==============================================================================
#This function will allow to read the given file and extract the data in liste
def Extract_data(file):
    listForLine=[] # liste for a ligne containing the element for the whole line 
    wholeListes=[] # liste containing all the smaller lists
    with open(file, 'r') as fil: 
        lines = fil.readlines()
        for line in lines:
            element= line.split() 
            #put the element of a line in a  list
            listForLine.append(element)

        for i in range (len(listForLine)):
            #select line that doesn't start with # or space to delet de comment of the gff/tsv file
            if ("#" and " " not in listForLine[i][0]):
                 wholeListes.append(listForLine[i])
        return  wholeListes 

def GeneDico(list):
    ListOfDicoGene=[]
    for element in list:
            dico_gene = {
            'chromosome':element[0],
            'source':element[1],
            'name':element[2],
            'start':element[3],
            'end':element[4],
            'score':element[5],
            'strand':element[6],
            'frame':element[7],
            'attribute':element[8]
            }
            ListOfDicoGene.append(dico_gene)
    print(ListOfDicoGene)
    return ListOfDicoGene

def TEDico(list):
    ListOfDicoTE=[]
    for element in list:
            dico_TE = {
            'chromosome':element[0],
            'source':element[1],
            'length_chromo':element[2],
            'name':element[3],
            'match':element[4],
            'start':element[5],
            'end':element[6],
            'Length':element[7],
            'Score':element[8],
            'Strand':element[9]
            }
            ListOfDicoTE.append(dico_TE)
    print(ListOfDicoTE)
    return ListOfDicoTE

def check_if_a_TE_in_gene(te,gene):#easiest version 
    return 

def check_if_a_TE_before_gene(te,gene):
    return

def check_if_a_TE_after_gene(te,gene):
    return


def writeDataOnFile(list):
    print("Please enter a file name .csv ")
    name=input()
    with open(name, 'w') as fil:
        print()
    fil.close()






#==============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==============================================================================

gene=Extract_data('Exemple_file/gene.tsv')
te=Extract_data('Exemple_file/te.tsv')

list_gene = GeneDico(gene)
list_te = TEDico(te)

print("liste gene ", list_gene)
print("liste te ", list_te)
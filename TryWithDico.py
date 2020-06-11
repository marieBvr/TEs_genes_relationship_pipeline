#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Ayse Ergun
# Script For TE project

#==============================================================================
#                           importations
#==============================================================================
import sys
import csv
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

        for i in range (1,len(listForLine)):
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
    #print(ListOfDicoGene)
    return ListOfDicoGene

def TEDico(list):
    ListOfDicoTE=[]
    for element in list:
            dico_TE = {
            'chromosome':element[0],
            'length_chromo':element[1],
            'name':element[2],
            'match':element[3],
            'start':element[4],
            'end':element[5],
            'length':element[6],
            'score':element[7],
            'strand':element[8]
            }
            ListOfDicoTE.append(dico_TE)
    #print(ListOfDicoTE)
    return ListOfDicoTE

def check_if_a_TE_in_gene(te,gene):
    for i in range(len(te)):
        for j in range(len(gene)):
            if(te[i]['start'] > gene[j]['start'] and te[i]['end']<gene[j]['end']):
                print(te[i]['name'] + ' is in the '+gene[j]['name']+' gene ! ')
    return 

def check_if_a_gene_in_TE(te,gene):
    return 

def check_if_a_TE_before_gene(te,gene):
    return

def check_if_a_TE_after_gene(te,gene):
    return

def calcul_distance():
    print('function calcul distance')


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

check_if_a_TE_in_gene(list_te, list_gene)

#test_push
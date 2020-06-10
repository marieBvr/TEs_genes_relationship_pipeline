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
            'seqname':element[0],
            'source':element[1],
            'feature':element[2],
            'start':element[3],
            'end':element[4],
            'score':element[5],
            'strand':element[6],
            'frame':element[7],
            'attribute':element[8]
            }
            ListOfDicoGene.append(dico_gene)
    return ListOfDicoGene

def TEDico(list):
    ListOfDicoTE=[]
    for element in list:
            dico_TE = {
            'seqname':element[0],
            'source':element[1],
            'feature':element[2],
            'start':element[3],
            'end':element[4],
            'score':element[5],
            'strand':element[6],
            'frame':element[7],
            'attribute':element[8]
            }
            ListOfDicoTE.append(dico_TE)
    return ListOfDicoTE






#==============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==============================================================================
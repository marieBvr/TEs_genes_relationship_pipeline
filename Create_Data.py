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
# https://www.w3schools.com/python/python_classes.asp
#==============================================================================
#                            functions
#==============================================================================
#This function will allow to read the given file and extract the data in liste
class Transposon :
    def __init__(self,startT, endT, nameT, strandT):
        self.startT=startT
        self.endT= endT
        self.nameT= nameT
        self.strandT=strandT

class Gene:
    def __init__(self,nameG, startG,endG,strandG):
        self.startG=startG
        self.endG= endG
        self.nameG= nameG
        self.strandG=strandG

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

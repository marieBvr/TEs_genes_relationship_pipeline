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

def Extract_data():



#==============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==============================================================================
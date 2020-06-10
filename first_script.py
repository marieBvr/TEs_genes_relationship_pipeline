#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#==============================================================================
#                           importations
#==============================================================================

import sys
import csv


#==============================================================================
#                            functions
#==============================================================================

def displayMenu():
    print("++++++++++++++++++++++++++++++++++")
    print("1- Show data from a gff file")
    print("2- Show only Start and End") #exemple of what we can do
    print("3- Put Start and End in a new file")
    print("4- Exit")
    print("++++++++++++++++++++++++++++++++++")

def readFile(file):
    fullList=[]
    listeliste=[]
    with open(file, 'r') as fil: 
        lines = fil.readlines()
        for line in lines:
            l= line.split() 
            #put the element of a line in a  list
            fullList.append(l)

        for i in range (len(fullList)):
            #select line that doesn't start with # to delet de comment of the gff file
            if ("#" and " " not in fullList[i][0]):
                listeliste.append(fullList[i])
        return listeliste 



def putDico(listeliste):
    ListOfDico=[]
    for liste in listeliste:
            dico_list = {
            'seqname':liste[0],
            'source':liste[1],
            'feature':liste[2],
            'start':liste[3],
            'end':liste[4],
            'score':liste[5],
            'strand':liste[6],
            'frame':liste[7],
            'attribute':liste[8]
            }
            ListOfDico.append(dico_list)
    return ListOfDico


def ShowSpecificData(ListOfDico):
    start=[]
    ends=[]
    returnListe=[]
    for i in range(len(ListOfDico)):
        start.append(ListOfDico[i]['start'])
        ends.append(ListOfDico[i]['end'])

    returnListe.append(["start"]) 
    returnListe.append(start)
    returnListe.append(["end"])
    returnListe.append(ends)   
    return returnListe


def writeOnFile(returnListe):
    print("Please enter a file name .csv ")
    name=input()
    with open(name, 'w') as fil:
        a=zip(returnListe[0],returnListe[2])
        b=zip(returnListe[1],returnListe[3])
        for i in returnListe:
            #much easier to put the list in a table 
            writer = csv.writer(fil, delimiter='\t')
            writer.writerows(a)
            writer.writerows(b)
        
    fil.close()
  





#==============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==============================================================================
print("-------------------------")
print("Welcome to this program" )
print("-------------------------")
print("What do you want to do ?")
print("-------------------------")

while(True):
    displayMenu()
    ans=int(input())
    if ans ==1 :
        listOflist=readFile('Marouch_3.1_braker214_PruarM.gff3')
        fullDico=putDico(listOflist)
        print(fullDico)

    elif ans ==2 :
        returnListe= ShowSpecificData(fullDico)
        print(returnListe)
    
    elif ans ==3 :
        writeOnFile(returnListe)

    elif ans==4 :
        sys.exit()





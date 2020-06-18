
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Ayse Ergun & Caroline Meguerditchian
# Script For TE project

#==============================================================================
#                           importations
#==============================================================================
import sys
import csv
import numpy as np
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
            'start': int(element[3]),
            'end': int(element[4]),
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
            'length_chromo':int(element[1]),
            'name':element[2],
            'match':element[3],
            'start': int(element[4]),
            'end': int(element[5]),
            'length':element[6],
            'score':element[7],
            'strand':element[8]
            }
            ListOfDicoTE.append(dico_TE)
    #print(ListOfDicoTE)
    return ListOfDicoTE

def check_superset_subset_genes(te,gene):
    for i in range(len(te)):
        #add a new column in the dictionnary 
        te[i]['superset_start'] = np.NAN
        te[i]['superset_end'] = np.NAN
        te[i]['subset_start'] = np.NAN
        te[i]['subset_end'] = np.NAN

        # loop to compare all the genes to the TE
        for j in range(len(gene)): 
            distances = calcul_distance(te[i],gene[j])

            #check if the TE is inside the gene
            if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] < 0):
                print(te[i]['name'], "is in", gene[j]['name'])
                te[i]['superset_start'] = gene[j]['start']
                te[i]['superset_end'] = gene[j]['end']
                te[i]['superset_end'] = gene[j]['end']
                
            
            #check if the TE is over the gene
            if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] > 0):
                print(te[i]['name'], "is over", gene[j]['name'])
                te[i]['subset_start'] = gene[j]['start']
                te[i]['subset_end'] = gene[j]['end']


def check_downstream_genes(te,gene):
    closest_gene = None

    for i in range(len(te)):
        te[i]['after_start'] = np.NAN
        te[i]['after_end'] = np.NAN
        te[i]['downstream_overlap'] = np.NAN


        for j in range(len(gene)): # loop to compare all the genes to the TE
            distances = calcul_distance(te[i],gene[j])

            #find downstream genes with overlap
            if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] > 0): 
                closest_gene = gene[j]
                te[i]['after_start'] = gene[j]['start']
                te[i]['after_end'] = gene[j]['end']
                te[i]['downstream_overlap'] = abs(distances[1])
                break

            #find genes downstream
            if(distances[0] < 0 and distances[1] > 0 and distances[2] > 0 and distances[3] > 0):
                closest_gene = gene[j]
                te[i]['after_start'] = gene[j]['start']
                te[i]['after_end'] = gene[j]['end']
                break

        #make sure that if the TE is followed by another TE there is no downstream gene
        for k in range(len(te)):
            start_value = te[k]['start']
            if(start_value > te[i]['end'] and start_value < closest_gene['start']):
                te[i]['after_start'] = np.NAN
                te[i]['after_end'] = np.NAN

def check_upstream_genes(te,gene):
    closest_gene = None
    gene = gene[::-1]

    for i in range(len(te)):
        te[i]['before_start'] = np.NAN
        te[i]['before_end'] = np.NAN
        te[i]['upstream_overlap'] = np.NAN

        for j in range(len(gene)): # loop to compare all the genes to the TE
            distances = calcul_distance(te[i],gene[j])
            
            #find overlap with gene upstream 
            if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0): 
                closest_gene = gene[j]
                te[i]['before_start'] = gene[j]['start']
                te[i]['before_end'] = gene[j]['end']
                te[i]['upstream_overlap'] = abs(distances[0])
                break

            #find genes upstream
            if(distances[0] > 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0):
                closest_gene = gene[j]
                te[i]['before_start'] = gene[j]['start']
                te[i]['before_end'] = gene[j]['end']
                break
    
        #make sure that if the TE is preceded by another TE there is no upstream gene
        for k in range(len(te)):
            if(te[i]['start'] > te[k]['end'] and closest_gene['end'] < te[k]['end']):
                te[i]['before_start'] = np.NAN
                te[i]['before_end'] = np.NAN

def calcul_distance(te,gene):
    distance1 = te['start'] - gene['end']
    distance2 = gene['start'] - te['end']
    distance3 = gene['end'] - te['end']
    distance4 = gene['start'] - te['start']
    
    return [distance1, distance2, distance3, distance4]


def writeDataOnFile(list_te):
    csv_content = []
    column_names = ["Start","End","before_start","before_end","after_start","after_end",
    "superset_start","superset_end","subset_start","subset_end","upstream_overlap","downstream_overlap"]
    for n in range(len(list_te)):
        csv_content.append([])
        csv_content[n].append(list_te[n]['start'])
        csv_content[n].append(list_te[n]['end'])
        csv_content[n].append(list_te[n]['before_start'])
        csv_content[n].append(list_te[n]['before_end'])
        csv_content[n].append(list_te[n]['after_start'])
        csv_content[n].append(list_te[n]['after_end'])
        csv_content[n].append(list_te[n]['superset_start'])
        csv_content[n].append(list_te[n]['superset_end'])
        csv_content[n].append(list_te[n]['subset_start'])
        csv_content[n].append(list_te[n]['subset_end'])
        csv_content[n].append(list_te[n]['upstream_overlap'])
        csv_content[n].append(list_te[n]['downstream_overlap'])
    with open('ResultFile.tsv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t')
        filewriter.writerow(column_names)
        for i in range(len(csv_content)):
            filewriter.writerow(csv_content[i])
    csvfile.close()






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

check_superset_subset_genes(list_te, list_gene)
check_downstream_genes(list_te, list_gene)
check_upstream_genes(list_te, list_gene)
writeDataOnFile(list_te)

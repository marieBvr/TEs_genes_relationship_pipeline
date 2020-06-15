
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
            'length_chromo':element[1],
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

            #check if the TE is inside the gene
            if(te[i]['start'] > gene[j]['start'] and te[i]['end']<gene[j]['end']):
                print(te[i]['name'], "is in", gene[j]['name'])
                te[i]['superset_start'] = gene[j]['start']
                te[i]['superset_end'] = gene[j]['end']
                te[i]['superset_end'] = gene[j]['end']
                
            
            #check if the TE is over the gene
            if (te[i]['start'] < gene[j]['start'] and te[i]['end'] > gene[j]['end']):
                print(te[i]['name'], "is over", gene[j]['name'])
                te[i]['subset_start'] = gene[j]['start']
                te[i]['subset_end'] = gene[j]['end']
    return 

def check_downstream_genes(te,gene):
    distance = 1000000 
    closest_gene = 0
    closest_gene_index = 0
    closest_gene_start = 0
    closest_gene_end = 0

    for i in range(len(te)):
        t_end = te[i]['end'] # End of the TE

        for j in range(len(gene)): # loop to compare all the genes to the TE
            g_end = gene[j]['end'] # End of the gene 
            
            #find genes that end after the end of the TE, and doesn't start before the start of the TE
            if(t_end < g_end and te[i]['start'] < gene[j]['start']): 

                #calculate the distance between gene and TE and stock the current gene if it is closer
                if g_end - t_end < distance: 
                    distance = g_end - t_end
                    #print(distance)
                    closest_gene_index = j #index of the new closest gene
                    closest_gene = gene[j]['name'] # id of the new closest gene
                    closest_gene_start = gene[j]['start'] #start of the new closest gene
                    closest_gene_end = gene[j]['end'] #end of the new closest gene

        #make sure that if the TE is followed by another TE there is no downstream gene
        for k in range(len(te)):
            start_value = te[k]['start']
            if(start_value > te[i]['end'] and start_value < gene[closest_gene_index]['start']):
                closest_gene = np.NAN
                closest_gene_start = np.NAN
                closest_gene_end = np.NAN

        distance = 100000

        te[i]['after_start'] = closest_gene_start
        te[i]['after_end'] = closest_gene_end
        print(closest_gene, "is downstream from ",te[i]['name'], closest_gene_start, closest_gene_end) # id du gène le plus proche du transposon d'intérêt

    return closest_gene_start, closest_gene_end

def check_upstream_genes(te,gene):
    distance = 1000000 
    closest_gene = 0
    closest_gene_index = 0
    closest_gene_start = 0
    closest_gene_end = 0

    for i in range(len(te)):
        t_start = te[i]['start'] # Start of the TE

        for j in range(len(gene)): # loop to compare all the genes to the TE
            g_start = gene[j]['start'] # start of the gene 
            
            #find genes that start before the start of the TE, and doesn't end after the end of the TE
            if(t_start > g_start and te[i]['end'] > gene[j]['end']): 

                #calculate the distance between gene and TE and stock the current gene if it is closer
                if t_start - g_start < distance: 
                    distance = t_start - g_start 
                    closest_gene_index = j #index of the new closest gene
                    closest_gene = gene[j]['name'] # id of the new closest gene
                    closest_gene_start = gene[j]['start'] #start of the new closest gene
                    closest_gene_end = gene[j]['end'] #end of the new closest gene
    
        #make sure that if the TE is preceded by another TE there is no upstream gene
        for k in range(len(te)):
            end_value = te[k]['end']
            if(end_value < te[i]['start'] and end_value > gene[closest_gene_index]['end']):
                closest_gene = np.NAN
                closest_gene_start = np.NAN
                closest_gene_end = np.NAN

        distance = 100000

        te[i]['before_start'] = closest_gene_start
        te[i]['before_end'] = closest_gene_end
        print(closest_gene, "is upstream from ",te[i]['name'], closest_gene_start, closest_gene_end) # id du gène le plus proche du transposon d'intérêt

    return closest_gene_start, closest_gene_end

def calcul_distance():
    print('function calcul distance')


def writeDataOnFile(list_te):
    csv_content = []
    column_names = ["Start","End","before_start","before_end","after_start","after_end",
    "superset_start","superset_end","subset_start","subset_end"]
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

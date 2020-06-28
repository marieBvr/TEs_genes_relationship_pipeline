
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
    i = 0
    chr_ = 'chr1'
    ListOfDicoGene=[[]]
    for element in list:
        if(chr_ != element[0]):
            chr_ = element[0]
            ListOfDicoGene.append([])
            i = i + 1
        dico_gene = {
        'chromosome':element[0],
        'source':element[1],
        'feature':element[2],
        'start': int(element[3]),
        'end': int(element[4]),
        'score':element[5],
        'strand':element[6],
        'frame':element[7],
        'attribute':element[8]
        }
        ListOfDicoGene[i].append(dico_gene)
    #print(ListOfDicoGene)
    return ListOfDicoGene

def TEDico(list):
    i = 0
    chr_ = 'chr1'
    ListOfDicoTE=[[]]
    for element in list:

        #make a new list for each chromosome
        if(chr_ != element[7]):
            chr_ = element[7]
            ListOfDicoTE.append([])
            i = i + 1
        dico_TE = {
        'type':element[12],
        'chromosome':element[7],
        'name':element[1],
        'start': int(element[8]),
        'end': int(element[9]),
        'score':element[15],
        'strand':element[10]
        }
        ListOfDicoTE[i].append(dico_TE)
    #print(ListOfDicoTE)
    return ListOfDicoTE

def check_superset_subset_genes(te,gene):
    #loop to look through each chromosome
    for ch in range(len(te)):

        #loop to look through each TE
        for i in range(len(te[ch])):
            #add a new column in the dictionnary 
            te[ch][i]['superset_feature'] = np.NAN
            te[ch][i]['superset_strand'] = np.NAN
            te[ch][i]['superset_start'] = np.NAN
            te[ch][i]['superset_end'] = np.NAN
            te[ch][i]['superset_id'] = np.NAN
            te[ch][i]['subset_strand'] = []
            te[ch][i]['subset_feature'] = []
            te[ch][i]['subset_start'] = []
            te[ch][i]['subset_end'] = []
            te[ch][i]['subset_id'] = []

            # loop to compare all the genes to the TE
            for j in range(len(gene[ch])): 
                distances = calcul_distance(te[ch][i],gene[ch][j])
                #print('this is j' , j)

                #check if the TE is inside the gene
                if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] < 0):
                    #print(te[i]['name'], "is in", gene[j]['name'])
                    te[ch][i]['superset_feature'] = gene[ch][j]['feature']
                    te[ch][i]['superset_strand'] = gene[ch][j]['strand']
                    te[ch][i]['superset_start'] = gene[ch][j]['start']
                    te[ch][i]['superset_end'] = gene[ch][j]['end']
                    te[ch][i]['superset_id'] = gene[ch][j]['attribute']
                    
                
                #check if the TE is over the gene
                if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] > 0):
                    #print(te[i]['name'], "is over", gene[j]['name'])
                    te[ch][i]['subset_strand'].append(gene[ch][j]['strand'])
                    te[ch][i]['subset_feature'].append(gene[ch][j]['feature'])
                    te[ch][i]['subset_start'].append(gene[ch][j]['start'])
                    te[ch][i]['subset_end'].append(gene[ch][j]['end'])
                    te[ch][i]['subset_id'].append(gene[ch][j]['attribute'])
            
            if(te[ch][i]['subset_start'] == []):
                te[ch][i]['subset_start'] = np.NAN
                te[ch][i]['subset_end'] = np.NAN
                te[ch][i]['subset_id'] = np.NAN
                te[ch][i]['subset_strand'] = np.NAN
                te[ch][i]['subset_feature'] = np.NAN


def check_downstream_genes(te,gene):
    closest_gene = None

    #loop to look through each chromosome
    for ch in range(len(te)):

        #loop to look through each TE
        for i in range(len(te[ch])):
            te[ch][i]['after_feature'] = np.NAN
            te[ch][i]['after_strand'] = np.NAN
            te[ch][i]['after_start'] = np.NAN
            te[ch][i]['after_end'] = np.NAN
            te[ch][i]['after_id'] = np.NAN
            te[ch][i]['downstream_overlap'] = np.NAN


            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calcul_distance(te[ch][i],gene[ch][j])
                te[ch][i]['Down_TEstart-Geneend'] = distances[0] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Down_Genestart-TEend'] = distances[1] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Down_Geneend-TEend'] = distances[2] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Down_Genestart-TEstart'] = distances[3] ######################################################## pour ajouter sur la distance sur le fichier

                #find downstream genes with overlap
                if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] > 0): 
                    closest_gene = gene[ch][j]
                    te[ch][i]['after_feature'] = gene[ch][j]['feature']
                    te[ch][i]['after_strand'] = gene[ch][j]['strand']
                    te[ch][i]['after_id'] = gene[ch][j]['attribute']
                    te[ch][i]['after_start'] = gene[ch][j]['start']
                    te[ch][i]['after_end'] = gene[ch][j]['end']
                    te[ch][i]['downstream_overlap'] = abs(distances[1])
                    break

                #find genes downstream
                if(distances[0] < 0 and distances[1] > 0 and distances[2] > 0 and distances[3] > 0):
                    closest_gene = gene[ch][j]
                    te[ch][i]['after_feature'] = gene[ch][j]['feature']
                    te[ch][i]['after_strand'] = gene[ch][j]['strand']
                    te[ch][i]['after_id'] = gene[ch][j]['attribute']
                    te[ch][i]['after_start'] = gene[ch][j]['start']
                    te[ch][i]['after_end'] = gene[ch][j]['end']
                    break

            #make sure that if the TE is followed by another TE there is no downstream gene
            for k in range(len(te[ch])):
                start_value = te[ch][k]['start']
                if(start_value > te[ch][i]['end'] and start_value < closest_gene['start']):
                    te[ch][i]['after_feature'] = np.NAN
                    te[ch][i]['after_strand'] = np.NAN
                    te[ch][i]['after_id'] = np.NAN
                    te[ch][i]['after_start'] = np.NAN
                    te[ch][i]['after_end'] = np.NAN

def check_upstream_genes(te,gene):
    closest_gene = None

    for ch in range(len(gene)): 
        gene[ch] = gene[ch][::-1]

    #loop to look through each chromosome
    for ch in range(len(te)):

        #loop to look through each TE
        for i in range(len(te[ch])):
            te[ch][i]['before_feature'] = np.NAN
            te[ch][i]['before_strand'] = np.NAN
            te[ch][i]['before_start'] = np.NAN
            te[ch][i]['before_end'] = np.NAN
            te[ch][i]['before_id'] = np.NAN
            te[ch][i]['upstream_overlap'] = np.NAN

            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calcul_distance(te[ch][i],gene[ch][j])
                te[ch][i]['Up_TEstart-Geneend'] = distances[0] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Up_Genestart-TEend'] = distances[1] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Up_Geneend-TEend'] = distances[2] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Up_Genestart-TEstart'] = distances[3] ######################################################## pour ajouter sur la distance sur le fichier

                #find overlap with gene upstream 
                if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0): 
                    closest_gene = gene[ch][j]
                    te[ch][i]['before_feature'] = gene[ch][j]['feature']
                    te[ch][i]['before_strand'] = gene[ch][j]['strand']
                    te[ch][i]['before_start'] = gene[ch][j]['start']
                    te[ch][i]['before_end'] = gene[ch][j]['end']
                    te[ch][i]['before_id'] = gene[ch][j]['attribute']
                    te[ch][i]['upstream_overlap'] = abs(distances[0])
                    break

                #find genes upstream
                if(distances[0] > 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0):
                    closest_gene = gene[ch][j]
                    te[ch][i]['before_feature'] = gene[ch][j]['feature']
                    te[ch][i]['before_strand'] = gene[ch][j]['strand']
                    te[ch][i]['before_start'] = gene[ch][j]['start']
                    te[ch][i]['before_end'] = gene[ch][j]['end']
                    te[ch][i]['before_id'] = gene[ch][j]['attribute']
                    break

            #make sure that if the TE is preceded by another TE there is no upstream gene
            for k in range(len(te[ch])):
                if(te[ch][i]['start'] > te[ch][k]['end'] and closest_gene['end'] < te[ch][k]['end']):
                    te[ch][i]['before_feature'] = np.NAN
                    te[ch][i]['before_strand'] = np.NAN
                    te[ch][i]['before_start'] = np.NAN
                    te[ch][i]['before_end'] = np.NAN
                    te[ch][i]['before_id'] = np.NAN

def calcul_distance(te,gene):
    distance1 = te['start'] - gene['end']
    distance2 = gene['start'] - te['end']
    distance3 = gene['end'] - te['end']
    distance4 = gene['start'] - te['start']
    
    return [distance1, distance2, distance3, distance4]


def writeDataOnFile(list_te):
    csv_content = []
    n = 0
    column_names = ["TE_Type","TE_id","chromosome","TE_strand","start","end","before_id",'before_feature',
    'before_strand','before_start','before_end',"after_id",'after_feature','after_strand',
    "after_start","after_end",'superset_id','superset_feature','superset_strand','superset_start',
    'superset_end','subset_id','subset_feature','subset_strand','subset_start','subset_end',
    'upstream_overlap',"downstream_overlap","Down_TEstart-Geneend","Down_Genestart-TEend",
    "Down_Geneend-TEend","Down_Genestart-TEstart",'Up_TEstart-Geneend','Up_Genestart-TEend',
    'Up_Geneend-TEend','Up_Genestart-TEstart']
    for c in range(len(list_te)):
        for t in range(len(list_te[c])):

            csv_content.append([])
            csv_content[n].append(list_te[c][t]['type'])
            csv_content[n].append(list_te[c][t]['name'])
            csv_content[n].append(list_te[c][t]['chromosome'])
            csv_content[n].append(list_te[c][t]['strand'])
            csv_content[n].append(list_te[c][t]['start'])
            csv_content[n].append(list_te[c][t]['end'])

            csv_content[n].append(list_te[c][t]['before_id'])
            csv_content[n].append(list_te[c][t]['before_feature'])
            csv_content[n].append(list_te[c][t]['before_strand'])
            csv_content[n].append(list_te[c][t]['before_start'])
            csv_content[n].append(list_te[c][t]['before_end'])


            csv_content[n].append(list_te[c][t]['after_id'])
            csv_content[n].append(list_te[c][t]['after_feature'])
            csv_content[n].append(list_te[c][t]['after_strand'])
            csv_content[n].append(list_te[c][t]['after_start'])
            csv_content[n].append(list_te[c][t]['after_end'])

            csv_content[n].append(list_te[c][t]['superset_id'])
            csv_content[n].append(list_te[c][t]['superset_feature'])
            csv_content[n].append(list_te[c][t]['superset_strand'])
            csv_content[n].append(list_te[c][t]['superset_start'])
            csv_content[n].append(list_te[c][t]['superset_end'])

            csv_content[n].append(list_te[c][t]['subset_id'])
            csv_content[n].append(list_te[c][t]['subset_feature'])
            csv_content[n].append(list_te[c][t]['subset_strand'])
            csv_content[n].append(list_te[c][t]['subset_start'])
            csv_content[n].append(list_te[c][t]['subset_end'])
            
            csv_content[n].append(list_te[c][t]['upstream_overlap'])
            csv_content[n].append(list_te[c][t]['downstream_overlap'])

            csv_content[n].append(list_te[c][t]['Down_TEstart-Geneend'])
            csv_content[n].append(list_te[c][t]['Down_Genestart-TEend'])
            csv_content[n].append(list_te[c][t]['Down_Geneend-TEend'])
            csv_content[n].append(list_te[c][t]['Down_Genestart-TEstart'])

            csv_content[n].append(list_te[c][t]['Up_TEstart-Geneend'])
            csv_content[n].append(list_te[c][t]['Up_Genestart-TEend'])
            csv_content[n].append(list_te[c][t]['Up_Geneend-TEend'])
            csv_content[n].append(list_te[c][t]['Up_Genestart-TEstart'])

            n = n + 1

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

gene=Extract_data(sys.argv[1])
te=Extract_data(sys.argv[2])
#gene=Extract_data('real_gene_data.tsv')
#te=Extract_data('real_LTR_data.tsv')
#print(te)

list_gene = GeneDico(gene)
list_te = TEDico(te)

#print("liste gene ", list_gene)
#print("liste te ", list_te)

check_superset_subset_genes(list_te, list_gene)
check_downstream_genes(list_te, list_gene)
check_upstream_genes(list_te, list_gene)
writeDataOnFile(list_te)
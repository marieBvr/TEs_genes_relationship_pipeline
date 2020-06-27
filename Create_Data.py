
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
            if ("#" or " " not in listForLine[i][0]):
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
    ListOfDicoTE=[]
    for element in list:
        #make a new list for each chromosome
        if(chr_ != element[0]):
            chr_ = element[0]
            ListOfDicoTE.append([])
            i = i + 1
        dico_TE = {
        'chr':element[0],
        'length_chromo':int(element[1]),
        'type':element[2],
        'match':element[3],
        'start': int(element[4]),
        'end': int(element[5]),
        'length':element[6],
        'score':element[7],
        'strand':element[8],
        'frame':element[9],
        'attribute':element[10],
        'code':element[11],
        'class':element[12],
        'TE_name':element[13],
        'TE_status':element[14]
        }
        ListOfDicoTE.append(dico_TE)
    #print(ListOfDicoTE)
    return ListOfDicoTE

def check_superset_subset_genes(te,gene):
    for ch in range(len(te)):

        #loop to look through each TE
        for i in range(len(te[ch])):
            #add a new column in the dictionnary 
            te[ch][i]['superset_feature'] = []
            te[ch][i]['superset_strand'] = []
            te[ch][i]['superset_start'] = []
            te[ch][i]['superset_end'] = []
            te[ch][i]['superset_id'] = []
            te[ch][i]['subset_strand'] = []
            te[ch][i]['subset_feature'] = []
            te[ch][i]['subset_start'] = []
            te[ch][i]['subset_end'] = []
            te[ch][i]['subset_id'] = []

        # loop to compare all the genes to the TE
            for j in range(len(gene[ch])): 
                distances = calcul_distance(te[ch][i],gene[ch][j])

                #check if the TE is inside the gene
                # ---------|    %%%%%%% gene %%%%%%%   |------------
                # ---------------|** TE ** | ------------------
                if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] < 0):
                    #print(te[ch][i]['name'], "is in", gene[ch][j]['name'])
                    te[ch][i]['superset_feature'].append(gene[ch][j]['feature'])##-> list
                    te[ch][i]['superset_strand'].append(gene[ch][j]['strand'])##-> list
                    te[ch][i]['superset_start'].append(gene[ch][j]['start'])##-> list
                    te[ch][i]['superset_end'].append(gene[ch][j]['end'])##-> list
                    te[ch][i]['superset_id'].append(gene[ch][j]['attribute'])##-> list
                    
                
                #check if the gene is inside the TE
                # ---------|    ****** TE *****     |------------
                # ---------------| %% gene %% | ------------------
                if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] > 0):
                    #print(te[ch][i]['name'], "is over", gene[ch][j]['name'])
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
    return te



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
            te[ch][i]['downstream_overlap_ID']=[]
            te[ch][i]['downstream_overlap_strand']=[]
            te[ch][i]['downstream_overlap_feature']=[]
            te[ch][i]['downstream_overlap_start']=[]
            te[ch][i]['downstream_overlap_end']=[]


            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calcul_distance(te[ch][i],gene[ch][j])
                te[ch][i]['Down_TEstart-Geneend'] = distances[0] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Down_Genestart-TEend'] = distances[1] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Down_Geneend-TEend'] = distances[2] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Down_Genestart-TEstart'] = distances[3] ######################################################## pour ajouter sur la distance sur le fichier

                #find downstream genes with overlap
                # ---------|    ****** TE *****     |------------
                # --------------------------- | %% gene %% | -----
                if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] > 0): 
                    closest_gene = gene[ch][j]
                    '''
                    te[ch][i]['after_feature'] = gene[ch][j]['feature']
                    te[ch][i]['after_strand'] = gene[ch][j]['strand']
                    te[ch][i]['after_id'] = gene[ch][j]['attribute']
                    te[ch][i]['after_start'] = gene[ch][j]['start']
                    te[ch][i]['after_end'] = gene[ch][j]['end']
                    '''
                    te[ch][i]['downstream_overlap'] = abs(distances[1])
                    te[ch][i]['downstream_overlap_ID'].append(gene[ch][j]['id'])
                    te[ch][i]['downstream_overlap_strand'].append(gene[ch][j]['strand'])
                    te[ch][i]['downstream_overlap_feature'].append(gene[ch][j]['feature'])
                    te[ch][i]['downstream_overlap_start'].append(gene[ch][j]['start'])
                    te[ch][i]['downstream_overlap_end'].append(gene[ch][j]['end'])
                    break

                #find genes downstream
                # ---------|    ****** TE *****     |----------------------------
                # -------------------------------------------- | %% gene %% | -----
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
    return te

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
            te[ch][i]['upstream_overlap_ID']=[]
            te[ch][i]['upstream_overlap_strand']=[]
            te[ch][i]['upstream_overlap_feature']=[]
            te[ch][i]['upstream_overlap_start']=[]
            te[ch][i]['upstream_overlap_end']=[]



            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calcul_distance(te[ch][i],gene[ch][j])
                te[ch][i]['Up_TEstart-Geneend'] = distances[0] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Up_Genestart-TEend'] = distances[1] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Up_Geneend-TEend'] = distances[2] ######################################################## pour ajouter sur la distance sur le fichier
                te[ch][i]['Up_Genestart-TEstart'] = distances[3] ######################################################## pour ajouter sur la distance sur le fichier

                #find overlap with gene upstream 
                # ---------|    ****** TE *****     |------------
                # --- | %% gene %% | ----------------------------
                if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0): 
                    closest_gene = gene[ch][j]
                    '''
                    te[ch][i]['before_feature'] = gene[ch][j]['feature']
                    te[ch][i]['before_strand'] = gene[ch][j]['strand']
                    te[ch][i]['before_start'] = gene[ch][j]['start']
                    te[ch][i]['before_end'] = gene[ch][j]['end']
                    te[ch][i]['before_id'] = gene[ch][j]['attribute']
                    '''
                    te[ch][i]['upstream_overlap'] = abs(distances[0])
                    
                    te[ch][i]['upstream_overlap_ID'].append(gene[ch][j]['id'])
                    te[ch][i]['upstream_overlap_strand'].append(gene[ch][j]['strand'])
                    te[ch][i]['upstream_overlap_feature'].append(gene[ch][j]['feature'])
                    te[ch][i]['upstream_overlap_start'].append(gene[ch][j]['start'])
                    te[ch][i]['upstream_overlap_end'].append(gene[ch][j]['end'])
                    break

                #find genes upstream
                # ----------------------------|    ****** TE *****     |------------
                # --- | %% gene %% | -----------------------------------------------
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
                    te[ch][i]['before_id'] = np.NaN

    print(te)
    return te

def calcul_distance(te,gene):
    distance1 = te['start'] - gene['end']
    distance2 = gene['start'] - te['end']
    distance3 = gene['end'] - te['end']
    distance4 = gene['start'] - te['start']
    
    return [distance1, distance2, distance3, distance4]


def writeDataOnFile(list_te):
    print(list_te)
    '''csv_content = []
    column_names = ["start","end","before_start","before_end","after_start","after_end",
    "superset_start","superset_end","subset_start","subset_end","upstream_overlap","downstream_overlap",
    "Down_TEstart-Geneend","Down_Genestart-TEend","Down_Geneend-TEend","Down_Genestart-TEstart",
    "Up_TEstart-Geneend","Up_Genestart-TEend","Up_Geneend-TEend","Up_Genestart-TEstart"]
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

        csv_content[n].append(list_te[n]['Down_TEstart-Geneend'])
        csv_content[n].append(list_te[n]['Down_Genestart-TEend'])
        csv_content[n].append(list_te[n]['Down_Geneend-TEend'])
        csv_content[n].append(list_te[n]['Down_Genestart-TEstart'])

        csv_content[n].append(list_te[n]['Up_TEstart-Geneend'])
        csv_content[n].append(list_te[n]['Up_Genestart-TEend'])
        csv_content[n].append(list_te[n]['Up_Geneend-TEend'])
        csv_content[n].append(list_te[n]['Up_Genestart-TEstart'])'''


    with open('ResultFile_TE.tsv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t')
        #filewriter.writerow(column_names)
        #for i in range(len(csv_content)):
        #    filewriter.writerow(csv_content[i])
        for i in range(len(list_te)):
            filewriter.writerow(list_te[i])
    csvfile.close()






#==============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==============================================================================

gene=Extract_data('real_gene_data.tsv')
te=Extract_data('real_TE_data.tsv')
print(te)
list_gene = GeneDico(gene)
list_te = TEDico(te)

#print("liste gene ", list_gene)
#print("liste te ", list_te)

liste_te1= check_superset_subset_genes(list_te, list_gene)
liste_te2=check_downstream_genes(liste_te1, list_gene)
liste_te3=check_upstream_genes(liste_te2, list_gene)
writeDataOnFile(liste_te3)

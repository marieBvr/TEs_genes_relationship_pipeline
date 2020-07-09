#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Ayse Ergun & Caroline Meguerditchian
# Script For TE project

#==============================================================================
#                               imports
#==============================================================================
import sys
import csv
import time
import numpy as np
import multiprocessing as mp
#==============================================================================
#                               functions
#==============================================================================
#This function will allow to read the given file and extract the data in list
def Extract_data(file):
    start_time = time.time()

    listForLine=[] # list for a line containing the element for the whole line 
    wholeLists=[] # list containing all the smaller lists
    with open(file, 'r') as fil: 
        lines = fil.readlines()
        for line in lines:
            element= line.split() 
            #put the element of a line in a  list
            listForLine.append(element)

        for i in range (1,len(listForLine)):
            #select line that doesn't start with # or space to delete the comment of the gff/tsv file
            if ("#" and " " not in listForLine[i][0]):
                 wholeLists.append(listForLine[i])

    elapsed_time = round((time.time() - start_time), 2)
    print("Extract_data time : ",elapsed_time)

    return  wholeLists 

def GeneDico(list):
    start_time = time.time()

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

    elapsed_time = round((time.time() - start_time), 2)
    print("GeneDico time : ",elapsed_time)

    return ListOfDicoGene

def TEDico(list):
    start_time = time.time()

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

    elapsed_time = round((time.time() - start_time), 2)
    print("TEDico time : ",elapsed_time)

    return ListOfDicoTE

def check_superset_subset_genes(queue,gene):
    start_time = time.time()

    te = queue.get()

    #loop to look through each chromosome
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
                distances = calculateDistance(te[ch][i],gene[ch][j])

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
                    
                
                #check if the TE is over the gene
                # ---------|    ****** TE *****     |------------
                # ---------------| %% gene %% | ------------------
                if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] > 0):
                    #print(te[i]['name'], "is over", gene[j]['name'])
                    te[ch][i]['subset_strand'].append(gene[ch][j]['strand'])
                    te[ch][i]['subset_feature'].append(gene[ch][j]['feature'])
                    te[ch][i]['subset_start'].append(gene[ch][j]['start'])
                    te[ch][i]['subset_end'].append(gene[ch][j]['end'])
                    te[ch][i]['subset_id'].append(gene[ch][j]['attribute'])
            
            #replace empty lists with NaN
            if(te[ch][i]['subset_start'] == []):
                te[ch][i]['subset_start'] = np.NAN
                te[ch][i]['subset_end'] = np.NAN
                te[ch][i]['subset_id'] = np.NAN
                te[ch][i]['subset_strand'] = np.NAN
                te[ch][i]['subset_feature'] = np.NAN
            if(te[ch][i]['superset_start'] == []):
                te[ch][i]['superset_feature'] = np.NAN
                te[ch][i]['superset_strand'] = np.NAN
                te[ch][i]['superset_start'] = np.NAN
                te[ch][i]['superset_end'] = np.NAN
                te[ch][i]['superset_id'] = np.NAN
    
    queue.put(te)

    elapsed_time = round((time.time() - start_time), 2)
    print("check_superset_subset_genes time : ",elapsed_time)


def check_downstream_genes(queue,gene):
    start_time = time.time()

    te = queue.get()

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

            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calculateDistance(te[ch][i],gene[ch][j])
                te[ch][i]['Down_TEstart-Geneend'] = distances[0] ######################################################## to add on the distance to the file
                te[ch][i]['Down_Genestart-TEend'] = distances[1] ######################################################## to add on the distance to the file
                te[ch][i]['Down_Geneend-TEend'] = distances[2] ######################################################## to add on the distance to the file
                te[ch][i]['Down_Genestart-TEstart'] = distances[3] ######################################################## to add on the distance to the file

                #find genes downstream
                # ---------|    ****** TE *****     |----------------------------
                # -------------------------------------------- | %% gene %% | -----
                if(distances[0] < 0 and distances[1] > 0 and distances[2] > 0 and distances[3] > 0 or distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] > 0):
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
                if(closest_gene != None):
                    if(start_value > te[ch][i]['end'] and start_value < closest_gene['start']):
                        te[ch][i]['after_feature'] = np.NAN
                        te[ch][i]['after_strand'] = np.NAN
                        te[ch][i]['after_id'] = np.NAN
                        te[ch][i]['after_start'] = np.NAN
                        te[ch][i]['after_end'] = np.NAN
                        te[ch][i]['Down_TEstart-Geneend'] = np.NAN
                        te[ch][i]['Down_Genestart-TEend'] = np.NAN
                        te[ch][i]['Down_Geneend-TEend'] = np.NAN
                        te[ch][i]['Down_Genestart-TEstart'] = np.NAN
    
    queue.put(te)

    elapsed_time = round((time.time() - start_time), 2)
    print("check_downstream_genes time : ",elapsed_time)

def check_upstream_genes(queue,gene):
    start_time = time.time()

    te = queue.get()

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

            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calculateDistance(te[ch][i],gene[ch][j])
                te[ch][i]['Up_TEstart-Geneend'] = distances[0] ######################################################## to add on the distance to the file
                te[ch][i]['Up_Genestart-TEend'] = distances[1] ######################################################## to add on the distance to the file
                te[ch][i]['Up_Geneend-TEend'] = distances[2] ######################################################## to add on the distance to the file
                te[ch][i]['Up_Genestart-TEstart'] = distances[3] ######################################################## to add on the distance to the file
                    

                #find genes upstream
                # ----------------------------|    ****** TE *****     |------------
                # --- | %% gene %% | -----------------------------------------------
                if(distances[0] > 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0 or distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0):
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
                    te[ch][i]['Up_TEstart-Geneend'] = np.NAN
                    te[ch][i]['Up_Genestart-TEend'] = np.NAN
                    te[ch][i]['Up_Geneend-TEend'] = np.NAN
                    te[ch][i]['Up_Genestart-TEstart'] = np.NAN

    queue.put(te)

    elapsed_time = round((time.time() - start_time), 2)
    print("check_upstream_genes time : ",elapsed_time)


def check_upstream_overlap(queue,gene):
    start_time = time.time()

    te = queue.get()

    #reverse the genes order
    for ch in range(len(gene)): 
        gene[ch] = gene[ch][::-1]

    #loop to look through each chromosome
    for ch in range(len(te)):

        #loop to look through each TE
        for i in range(len(te[ch])):
            te[ch][i]['upstream_overlap'] = []
            te[ch][i]['upstream_overlap_ID']=[]
            te[ch][i]['upstream_overlap_strand']=[]
            te[ch][i]['upstream_overlap_feature']=[]
            te[ch][i]['upstream_overlap_start']=[]
            te[ch][i]['upstream_overlap_end']=[]

            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calculateDistance(te[ch][i],gene[ch][j])

            #find overlap with gene upstream 
                # ---------|    ****** TE *****     |------------
                # --- | %% gene %% | ----------------------------
                if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0): 
                    te[ch][i]['upstream_overlap'].append(abs(distances[0]))                 
                    te[ch][i]['upstream_overlap_ID'].append(gene[ch][j]['attribute'])
                    te[ch][i]['upstream_overlap_strand'].append(gene[ch][j]['strand'])
                    te[ch][i]['upstream_overlap_feature'].append(gene[ch][j]['feature'])
                    te[ch][i]['upstream_overlap_start'].append(gene[ch][j]['start'])
                    te[ch][i]['upstream_overlap_end'].append(gene[ch][j]['end'])

    queue.put(te)

    elapsed_time = round((time.time() - start_time), 2)
    print("check_upstream_overlap time : ",elapsed_time)

def check_downstream_overlap(queue,gene):
    start_time = time.time()

    te = queue.get()

    #loop to look through each chromosome
    for ch in range(len(te)):

        #loop to look through each TE
        for i in range(len(te[ch])):
            te[ch][i]['downstream_overlap'] = []
            te[ch][i]['downstream_overlap_ID']=[]
            te[ch][i]['downstream_overlap_strand']=[]
            te[ch][i]['downstream_overlap_feature']=[]
            te[ch][i]['downstream_overlap_start']=[]
            te[ch][i]['downstream_overlap_end']=[]


            for j in range(len(gene[ch])): # loop to compare all the genes to the TE
                distances = calculateDistance(te[ch][i],gene[ch][j])

                #find downstream genes with overlap
                # ---------|    ****** TE *****     |------------
                # --------------------------- | %% gene %% | -----
                if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] > 0):                     
                    te[ch][i]['downstream_overlap'].append(abs(distances[1]))
                    te[ch][i]['downstream_overlap_ID'].append(gene[ch][j]['attribute'])
                    te[ch][i]['downstream_overlap_strand'].append(gene[ch][j]['strand'])
                    te[ch][i]['downstream_overlap_feature'].append(gene[ch][j]['feature'])
                    te[ch][i]['downstream_overlap_start'].append(gene[ch][j]['start'])
                    te[ch][i]['downstream_overlap_end'].append(gene[ch][j]['end'])

    queue.put(te)

    elapsed_time = round((time.time() - start_time), 2)
    print("check_downstream_overlap time : ",elapsed_time)
    


def calculateDistance(te,gene):
    distance1 = te['start'] - gene['end']
    distance2 = gene['start'] - te['end']
    distance3 = gene['end'] - te['end']
    distance4 = gene['start'] - te['start']
    
    return [distance1, distance2, distance3, distance4]


def writeDataOnFile(list_te):
    start_time = time.time()

    csv_content = []
    n = 0
    column_names = ["TE_Type","TE_id","chromosome","TE_strand","start","end","before_id",'before_feature',
    'before_strand','before_start','before_end',"after_id",'after_feature','after_strand',
    "after_start","after_end",'superset_id','superset_feature','superset_strand','superset_start',
    'superset_end','subset_id','subset_feature','subset_strand','subset_start','subset_end',
    'upstream_overlap','upstream_overlap_ID','upstream_overlap_strand','upstream_overlap_feature',
    'upstream_overlap_start','upstream_overlap_end',"downstream_overlap",
    'downstream_overlap_ID','downstream_overlap_strand','downstream_overlap_feature',
    'downstream_overlap_start','downstream_overlap_end',"Down_TEstart-Geneend","Down_Genestart-TEend",
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
            csv_content[n].append(list_te[c][t]['upstream_overlap_ID'])
            csv_content[n].append(list_te[c][t]['upstream_overlap_strand'])
            csv_content[n].append(list_te[c][t]['upstream_overlap_feature'])
            csv_content[n].append(list_te[c][t]['upstream_overlap_start'])
            csv_content[n].append(list_te[c][t]['upstream_overlap_end'])

            csv_content[n].append(list_te[c][t]['downstream_overlap'])
            csv_content[n].append(list_te[c][t]['downstream_overlap_ID'])
            csv_content[n].append(list_te[c][t]['downstream_overlap_strand'])
            csv_content[n].append(list_te[c][t]['downstream_overlap_feature'])
            csv_content[n].append(list_te[c][t]['downstream_overlap_start'])
            csv_content[n].append(list_te[c][t]['downstream_overlap_end'])

            csv_content[n].append(list_te[c][t]['Down_TEstart-Geneend'])
            csv_content[n].append(list_te[c][t]['Down_Genestart-TEend'])
            csv_content[n].append(list_te[c][t]['Down_Geneend-TEend'])
            csv_content[n].append(list_te[c][t]['Down_Genestart-TEstart'])

            csv_content[n].append(list_te[c][t]['Up_TEstart-Geneend'])
            csv_content[n].append(list_te[c][t]['Up_Genestart-TEend'])
            csv_content[n].append(list_te[c][t]['Up_Geneend-TEend'])
            csv_content[n].append(list_te[c][t]['Up_Genestart-TEstart'])

            n = n + 1

    with open('ResultFile_LTR_multi.tsv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t')
        filewriter.writerow(column_names)
        for i in range(len(csv_content)):
            filewriter.writerow(csv_content[i])
    csvfile.close()

    elapsed_time = round((time.time() - start_time), 2)
    print("writeDataOnFile time : ",elapsed_time)
    



#==============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==============================================================================

gene=Extract_data('real_gene_data.tsv')
te=Extract_data('real_TE_data.tsv')

#print(te)
list_gene = GeneDico(gene)
list_te = TEDico(te)

#print("list gene ", list_gene)
#print("list te ", list_te)


if __name__ == '__main__':

    #creating the queues to give the data to the functions and get the returned data
    queue1 = mp.Queue()
    queue2 = mp.Queue()
    queue3 = mp.Queue()
    queue4 = mp.Queue()
    queue5 = mp.Queue()

    #creating the process for each function
    a = mp.Process(target=check_superset_subset_genes,args=(queue1,list_gene))
    b = mp.Process(target=check_downstream_genes,args=(queue2,list_gene))
    c = mp.Process(target=check_upstream_genes,args=(queue3,list_gene))
    d = mp.Process(target=check_upstream_overlap,args=(queue4,list_gene))
    e = mp.Process(target=check_downstream_overlap,args=(queue5,list_gene))
    
    #starting all process
    a.start()
    b.start()
    c.start()
    d.start()
    e.start()

    #putting the list of TE in the queue for the process to use them
    queue1.put(list_te)
    queue2.put(list_te)
    queue3.put(list_te)
    queue4.put(list_te)
    queue5.put(list_te)

    #getting the data back from each function
    l1 = queue1.get()
    l2 = queue2.get()
    l3 = queue3.get()
    l4 = queue4.get()
    l5 = queue5.get()

    #finishing every process
    a.join()
    b.join()
    c.join()
    d.join()
    e.join()

    lists = [l2,l3,l4,l5]

    #loop to reconstitute the final list from every list given by each function
    for l in range(len(lists)):
        for i in range(len(l1)):
            length = len(lists[l][i])
            for j in range(length):
                keys = list(lists[l][i][j].keys())
                for key in keys:
                    l1[i][j][key] = lists[l][i][j][key]

    writeDataOnFile(l1)

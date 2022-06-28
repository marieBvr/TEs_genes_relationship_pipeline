#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Ayse Ergun & Caroline Meguerditchian
# Script for TE project

#==============================================================================
#                           importations
#==============================================================================
import sys
import csv
import numpy as np
import time
import multiprocessing as mp
import argparse
#==============================================================================
#                            functions
#==============================================================================
def Extract_gene_data(file):
    """
    This function reads the given file and extract the data in list
    """
    start_time = time.time()
    header = 0
    i = 0
    genes_list = []
    with open(file, 'r') as fil: 
        lines = fil.readlines()
        for line in lines:
            if header == 0:
                header += 1
                continue
            else:
                element = line.split("\t")
                if element[2] == "gene":
                    item_gene = {
                        'chr':element[0],
                        'source':element[1],
                        'feature':element[2],
                        'start': int(element[3]),
                        'end': int(element[4]),
                        'score':element[5],
                        'strand':element[6],
                        'phase':element[7],
                        'id':element[8],
                        'attribute':element[9]
                        }
                    genes_list.append(item_gene)
                    i += 1

    elapsed_time = round((time.time() - start_time), 2)
    print("Extract_data time : ", elapsed_time)

    return  genes_list


def Extract_TEs_data(file):
    """
    Generate a list of dict with specific column names
    """
    start_time = time.time()
    header = 0
    i = 0
    TEs_list = []
    with open(file, 'r') as fil: 
        lines = fil.readlines()
        for line in lines:
            if header == 0:
                header += 1
                continue
            else:
                element = line.split('\t')
                item_TE = {
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
                TEs_list.append(item_TE)
                i += 1

    elapsed_time = round((time.time() - start_time), 2)
    print("Extract_data time : ", elapsed_time)

    return  TEs_list


def check_downstream(genes_list, TEs_list):
    """
    Look for the downstream TE for each gene
    """
    start_time = time.time()
    result = []
    item = {}
    for gene in genes_list:
        closest_TE = None
        for TE in TEs_list:
            if TE['chr'] == gene['chr']:
                distances = calcul_distance(TE, gene)
                if(distances[0] > 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0):
                    item = {
                        "Chr": TE['chr'],
                        "TE_type": TE['type'],
                        "TE_id": TE['attribute'],
                        "TE_strand": TE['strand'],
                        "TE_code": TE['code'],
                        "TE_start": TE['start'],
                        "TE_end": TE['end'],
                        "Gene_id": gene['id'],
                        "Gene_strand": gene['strand'],
                        "Gene_start": gene['start'],
                        "Gene_end": gene['end'],
                        "TEstart - GeneEnd": distances[0],
                        "GeneStart - TEend": distances[1],
                        "GeneEnd - TEend": distances[2],
                        "GeneStart - TEstart": distances[3],
                        "relationship": 'downstream'
                    
                    }
                    if closest_TE == None:
                        closest_TE = item
                    elif int(distances[0]) < int(closest_TE["TEstart - GeneEnd"]):
                        closest_TE = item
        if closest_TE != None:
            result.append(closest_TE)
    elapsed_time = round((time.time() - start_time), 2)
    print("check_downstream_genes time : ",elapsed_time)
    return result

def _check_upstream(genes_list, TEs_list):
    """
    Look for the upstream TE for each Gene
    """
    start_time = time.time()
    result = []
    item = {}
    for gene in genes_list:
        closest_TE = None
        for TE in TEs_list:
            if TE['chr'] == gene['chr']:
                distances = calcul_distance(TE, gene)
                if(distances[0] < 0 and distances[1] > 0 and distances[2] > 0 and distances[3] > 0):
                    item = {
                        "Chr": TE['chr'],
                        "TE_type": TE['type'],
                        "TE_id": TE['attribute'],
                        "TE_strand": TE['strand'],
                        "TE_code": TE['code'],
                        "TE_start": TE['start'],
                        "TE_end": TE['end'],
                        "Gene_id": gene['id'],
                        "Gene_strand": gene['strand'],
                        "Gene_start": gene['start'],
                        "Gene_end": gene['end'],
                        "TEstart - GeneEnd": distances[0],
                        "GeneStart - TEend": distances[1],
                        "GeneEnd - TEend": distances[2],
                        "GeneStart - TEstart": distances[3],
                        "relationship": 'upstream'
                    
                    }
                    if closest_TE == None:
                        closest_TE = item
                    elif int(distances[1]) < int(closest_TE["GeneStart - TEend"]):
                        closest_TE = item
        if closest_TE != None:
            result.append(closest_TE)
    elapsed_time = round((time.time() - start_time), 2)
    print("check_upstream_genes time : ", elapsed_time)
    return result
    
    
def _check_upstream_overlap(genes_list, TEs_list):
    """
    Look for the upstream TE for each Gene
    """
    start_time = time.time()
    result = []
    item = {}
    to_remove = []
    for gene in genes_list:
        for TE in TEs_list:
            if TE['chr'] == gene['chr']:
                distances = calcul_distance(TE, gene)
                if(distances[0] < 0 and distances[1] < 0 and distances[2] > 0 and distances[3] > 0): 
                    item = {
                        "Chr": TE['chr'],
                        "TE_type": TE['type'],
                        "TE_id": TE['attribute'],
                        "TE_strand": TE['strand'],
                        "TE_code": TE['code'],
                        "TE_start": TE['start'],
                        "TE_end": TE['end'],
                        "Gene_id": gene['id'],
                        "Gene_strand": gene['strand'],
                        "Gene_start": gene['start'],
                        "Gene_end": gene['end'],
                        "TEstart - GeneEnd": distances[0],
                        "GeneStart - TEend": distances[1],
                        "GeneEnd - TEend": distances[2],
                        "GeneStart - TEstart": distances[3],
                        "relationship": 'upstream_overlap'
                    
                    }
                    result.append(item)
                    to_remove.append(gene)
    for r in to_remove:
        if r in genes_list:
            genes_list.remove(r)
    elapsed_time = round((time.time() - start_time), 2)
    print("check_upstream_genes time : ", elapsed_time)
    return [result, genes_list]



def check_downstream_overlap(genes_list, TEs_list):
    """
    Look for the upstream TE for each Gene
    """
    start_time = time.time()
    result = []
    item = {}
    to_remove = []
    for gene in genes_list:
        for TE in TEs_list:
            if TE['chr'] == gene['chr']:
                distances = calcul_distance(TE, gene)
                if(distances[0] < 0 and distances[1] < 0 and distances[2] < 0 and distances[3] < 0):  
                    item = {
                        "Chr": TE['chr'],
                        "TE_type": TE['type'],
                        "TE_id": TE['attribute'],
                        "TE_strand": TE['strand'],
                        "TE_code": TE['code'],
                        "TE_start": TE['start'],
                        "TE_end": TE['end'],
                        "Gene_id": gene['id'],
                        "Gene_strand": gene['strand'],
                        "Gene_start": gene['start'],
                        "Gene_end": gene['end'],
                        "TEstart - GeneEnd": distances[0],
                        "GeneStart - TEend": distances[1],
                        "GeneEnd - TEend": distances[2],
                        "GeneStart - TEstart": distances[3],
                        "relationship": 'downstream_overlap'
                    
                    }
                    result.append(item)
                    to_remove.append(gene)
    for r in to_remove:
        if r in genes_list:
            genes_list.remove(r)
    elapsed_time = round((time.time() - start_time), 2)
    print("check_upstream_genes time : ", elapsed_time)
    return [result, genes_list]
    
    
def check_subset_superset(genes_list, TEs_list):
    """
    Look for the upstream TE for each Gene
    """
    start_time = time.time()
    result = []
    item = {}
    for gene in genes_list:
        for TE in TEs_list:
            if TE['chr'] == gene['chr']:
                distances = calcul_distance(TE, gene)
                item = None
                if(distances[0] < 0 and distances[1] < 0 and distances[2] >= 0 and distances[3] <= 0): # subset 
                    item = {
                        "Chr": TE['chr'],
                        "TE_type": TE['type'],
                        "TE_id": TE['attribute'],
                        "TE_strand": TE['strand'],
                        "TE_code": TE['code'],
                        "TE_start": TE['start'],
                        "TE_end": TE['end'],
                        "Gene_id": gene['id'],
                        "Gene_strand": gene['strand'],
                        "Gene_start": gene['start'],
                        "Gene_end": gene['end'],
                        "TEstart - GeneEnd": distances[0],
                        "GeneStart - TEend": distances[1],
                        "GeneEnd - TEend": distances[2],
                        "GeneStart - TEstart": distances[3],
                        "relationship": 'subset'
                    
                    }
                elif(distances[0] < 0 and distances[1] < 0 and distances[2] <= 0 and distances[3] >= 0): #superset 
                    item = {
                        "Chr": TE['chr'],
                        "TE_type": TE['type'],
                        "TE_id": TE['attribute'],
                        "TE_strand": TE['strand'],
                        "TE_code": TE['code'],
                        "TE_start": TE['start'],
                        "TE_end": TE['end'],
                        "Gene_id": gene['id'],
                        "Gene_strand": gene['strand'],
                        "Gene_start": gene['start'],
                        "Gene_end": gene['end'],
                        "TEstart - GeneEnd": distances[0],
                        "GeneStart - TEend": distances[1],
                        "GeneEnd - TEend": distances[2],
                        "GeneStart - TEstart": distances[3],
                        "relationship": 'superset'
                    
                    }
                if item != None:
                    result.append(item)
    elapsed_time = round((time.time() - start_time), 2)
    print("check_upstream_genes time : ", elapsed_time)
    return result


def calcul_distance(te, gene):
    distance1 = te['start'] - gene['end']
    distance2 = gene['start'] - te['end']
    distance3 = gene['end'] - te['end']
    distance4 = gene['start'] - te['start']
    return [distance1, distance2, distance3, distance4]
    
    
def write_results(final_results, output_file):
    """
    Write the output file
    """
    start_time=time.time()
    column_names = ['Chr', 'Type', 'TE_id','TE_strand','TE_code', 'TE_start','TE_end','Gene_id',
    'Gene_strand','Gene_start','Gene_end','TEstart-Geneend','Genestart-TEend',
    'Geneend-TEend','Genestart-TEstart','relationship']
    with open(output_file, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t')
        filewriter.writerow(column_names)
        for item in final_results:
            # print(item.values())
            filewriter.writerow(item.values())
    csvfile.close()
    elapsed_time = round((time.time() - start_time), 2)
    print("writeDataOnFile time : ",elapsed_time)


def _set_options():
    """
    Define program options
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene',help='Genes position file as GFF.', action='store',required=True, type=str,dest='geneFile')
    parser.add_argument('-te', '--transposable_element', help='TEs position file as GFF.', action='store',required=True, type=str,dest='teFile')
    parser.add_argument('-o', '--out', help='The output file.', action='store', type=str, default='Resting_result.tsv', dest='output')
    args = parser.parse_args()
    return args

#==============================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#==============================================================================
if __name__ == '__main__':
    # Retrieve program options
    options = _set_options()
    genes = Extract_gene_data(options.geneFile)
    TEs = Extract_TEs_data(options.teFile)
    back_genes = genes.copy()

    # downstream
    results_d = check_downstream_overlap(genes, TEs)
    downstream_overlap = results_d[0]
    genes = results_d[1]
    downstream  = check_downstream(genes, TEs)

    # upstream
    results_up = _check_upstream_overlap(back_genes, TEs)
    upstream_overlap = results_up[0]
    genes = results_up[1]
    upstream  = _check_upstream(genes, TEs)

    subset_superset = check_subset_superset(back_genes, TEs)
    final_results = np.concatenate((
        downstream,
        upstream,
        downstream_overlap,
        upstream_overlap,
        subset_superset
        ))
    write_results(final_results, options.output)

    #creating the queues to give the data to the functions and get the returned data
    # queue1 = mp.Queue()
    # queue2 = mp.Queue()
    # queue3 = mp.Queue()
    # queue4 = mp.Queue()
    # queue5 = mp.Queue()
	
    #creating the process for each function
    # a = mp.Process(target=check_superset_subset_genes,args=(queue1,list_gene))
    # b = mp.Process(target=check_downstream_genes,args=(queue2,list_gene))
    # c = mp.Process(target=check_upstream_genes,args=(queue3,list_gene))
    # d = mp.Process(target=check_upstream_overlap,args=(queue4,list_gene))
    # e = mp.Process(target=check_downstream_overlap,args=(queue5,list_gene))
	
    #starting all process
    # a.start()
    # b.start()
    # c.start()
    # d.start()
    # e.start()
	
    #putting the list of TE in the queue for the process to use them
    # queue1.put(list_te)
    # queue2.put(list_te)
    # queue3.put(list_te)
    # queue4.put(list_te)
    # queue5.put(list_te)
	
    #getting the data back from each function
    # l1 = queue1.get()
    # l2 = queue2.get()
    # l3 = queue3.get()
    # l4 = queue4.get()
    # l5 = queue5.get()
	
    #finishing every process
    # a.join()
    # b.join()
    # c.join()
    # d.join()
    # e.join()
    # listes = [l2,l3,l4,l5]
    #loop to reconstitute the final list from every list given by each function
    # for l in range(len(listes)):
        # for i in range(len(l1)):
            # length = len(listes[l][i])
            # for j in range(length):
                # keys = list(listes[l][i][j].keys())
                # for key in keys:
                    # l1[i][j][key] = listes[l][i][j][key]
    # writeDataOnFile(l1)
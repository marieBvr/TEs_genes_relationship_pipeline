#!/usr/bin/env python3

# After running:
# python3 Multiprocessing/Create_Data_multipro_reformatted.py \
# -g data/Gene_testing_data.tsv \
# -te data/Transposon_testing_data.tsv \
# -o result/output_TE.tsv

# This script aims at taking into account the strand of genes
# in case gene on strand (-), TEs downstream of the gene
# become up and vice versa.

import argparse
import csv

dict_relationship = {
    "downstream": "upstream",
    "downstream_overlap": "upstream_overlap",
    "upstream": "downstream",
    "upstream_overlap": "downstream_overlap"
}

def strand_correction(inputFile, outputFile):
    header = 0
    with open(inputFile, 'r') as input:
        lines = input.readlines()
        with open(outputFile, 'w') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t')
            for line in lines:
                if header == 0:
                    header += 1
                    # write header
                    filewriter.writerow(line.strip('\n').split("\t"))
                    continue
                else:
                    # remove \n at the end of the line
                    element = line.strip('\n').split("\t")
                    gene_strand = element[8]
                    # if (-) change relationship
                    if gene_strand == "-":
                        relation = element[len(element)-1].replace("\n", "")
                        element[len(element)-1] = dict_relationship[relation]
                        filewriter.writerow(element)
                    else:
                        filewriter.writerow(element)
        csvfile.close()


def _set_options():
    """
    Define program options
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',help='Output of previous TE_multipro_refromatted script.', action='store', required=True, type=str,dest='input')
    parser.add_argument('-o', '--out', help='The output file.', action='store', type=str, default='post_result.tsv', dest='output')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # Retrieve program options
    options = _set_options()
    strand_correction(options.input, options.output)

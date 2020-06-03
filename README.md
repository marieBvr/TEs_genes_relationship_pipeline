# TE_abricot

You can see the report progression with the following link : https://fr.overleaf.com/read/qwygcwcjmxdd  

Distribution of TEs and their relationship to genes (and gene expression) in Apricot

Calculate the relationship between TEs and genes in basing on the TE annotation gff (REPET) and genome annotation gff (Braker214):



1)	TEs overlap with Genes
2)	Distance between TE with upstream Gene and downstream Gene

To do 
1.	Cross the interested information between TE gff (file: info_reswagMa3S_refTEsgff.tsv ) and Gene gff (file: Marouch_3.1_braker214_PruarM.gff3) to create a big table by: python scripts and/or Bedtools closest or Bedtools intersect.
Advise: output = keep all columns from TE gff file and add columns from Gene gff (chr	transdecoder	mRNA	start	end	.	+(strand)	0), and add new column(s) for the distance calculation between TE and Gene (The details of sense or antisense can be generated in the second step with R from this output)
2.	Generate the interested information in details (all questions above) with R (visualize by graphs, always “write table” of each graph). 

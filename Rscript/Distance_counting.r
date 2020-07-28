library(ggplot2)

result_file = read.csv(file = 'ResultFile_LTR_Mands.csv', sep = '\t', header = TRUE)


#To convert python list into vector
#library(stringr)
#a = str_remove_all(result_file2$subset_id[1], "[\"',\\[\\]]")
#b = strsplit(a, "\\s+")[[1]]


count_gene_upstream = function(result,sens,distance_min,distance_max){
count = 0
for(i in 1:length(result$TE_id)){
  upstream_gene_end = result$before_end[i]
  if(is.nan(upstream_gene_end) == FALSE & result$TE_strand[i] == sens & result$Up_TEstart.Geneend[i] > distance_min & result$Up_TEstart.Geneend[i] < distance_max){
    count = count + 1
    }
  }
  #print(count)
  return(count)
}

count_gene_downstream = function(result,sens,distance_min,distance_max){
  count = 0
  for(i in 1:length(result$TE_id)){
    downstream_gene_end = result$after_end[i]
    if(is.nan(downstream_gene_end) == FALSE & result$TE_strand[i] == sens & result$Down_Genestart.TEend[i] > distance_min & result$Down_Genestart.TEend[i] < distance_max){
      count = count + 1
    }
  }
  #print(count)
  return(count)
}

c1 = count_gene_upstream(result_file,'-',0,500)
c2 = count_gene_upstream(result_file,'-',501,1000)
c3 = count_gene_upstream(result_file,'-',1001,2000)
c4 = count_gene_upstream(result_file,'-',2001,max(result_file$end))

c5 = count_gene_upstream(result_file,'+',0,500)
c6 = count_gene_upstream(result_file,'+',501,1000)
c7 = count_gene_upstream(result_file,'+',1001,2000)
c8 = count_gene_upstream(result_file,'+',2001,max(result_file$end))

upstream_val = c(c5,c1,c6,c2,c7,c3,c8,c4)

c9 = count_gene_downstream(result_file,'-',0,500)
c10 = count_gene_downstream(result_file,'-',501,1000)
c11 = count_gene_downstream(result_file,'-',1001,2000)
c12 = count_gene_downstream(result_file,'-',2001,max(result_file$end))

c13 = count_gene_downstream(result_file,'+',0,500)
c14 = count_gene_downstream(result_file,'+',501,1000)
c15 = count_gene_downstream(result_file,'+',1001,2000)
c16 = count_gene_downstream(result_file,'+',2001,max(result_file$end))

downstream_val = c(c13,c9,c14,c10,c15,c11,c16,c12)

#graphique ggplot upstream
distance = c(rep("0-0.5kbp" , 2), rep("0.5-1kbp" , 2), rep("1-2kbp", 2), rep("2kbp+", 2))
strand = rep(c("+","-"),4)
value = upstream_val
data = data.frame(distance,strand,value)
plot1 = ggplot(data, aes(fill=strand, y=value, x=distance )) + geom_bar(position="dodge", stat="identity") + ggtitle("Number of upstream genes adjacent to LTR") + ylab("number of gene") + xlab("distance to the LTR")
print(plot1)


#graphique ggplot downstream
distance = c(rep("0-0.5kbp" , 2), rep("0.5-1kbp" , 2), rep("1-2kbp", 2), rep("2kbp+", 2))
strand = rep(c("+","-"),4)
value = downstream_val
data = data.frame(distance,strand,value)
plot2 = ggplot(data, aes(fill=strand, y=value, x=distance )) + geom_bar(position="dodge", stat="identity") + ggtitle("Number of downstream genes adjacent to LTR") + ylab("number of gene") + xlab("distance to the LTR")
print(plot2)

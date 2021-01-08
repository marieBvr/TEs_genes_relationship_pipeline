if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
}
if(!require("optparse")){
    install.packages("optparse")
    library(optparse)
}

# script options
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name i.e. 'Resting_result/output_LTR.tsv'", metavar="character"),
  make_option(c("-p", "--pdf"), type="character", default="Resting_result/distance_TE_results.pdf",
              help="output filename name (PDF) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Resting_result/distance_TE_results.csv",
              help="output filename name (CSV) [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# no option given to script
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file). ", call.=FALSE)
}

#load result file
result_file = read.csv(file = opt$file, sep = '\t', header = TRUE)

#count genes in the before gene column
count_gene_upstream = function(result,distance_min,distance_max){
  count_up_plus_strand = 0
  count_down_minus_strand = 0
  for(i in 1:length(result$TE_id)){
    upstream_gene_end = result$before_end[i]
    if(is.nan(upstream_gene_end) == FALSE & result$TE_strand[i] == '+' & result$Up_TEstart.Geneend[i] > distance_min & result$Up_TEstart.Geneend[i] < distance_max){
      count_up_plus_strand = count_up_plus_strand + 1
    }
    if(is.nan(upstream_gene_end) == FALSE & result$TE_strand[i] == '-' & result$Up_TEstart.Geneend[i] > distance_min & result$Up_TEstart.Geneend[i] < distance_max){
      count_down_minus_strand = count_down_minus_strand + 1
    }
  }
  count = c(count_up_plus_strand,count_down_minus_strand)
  return(count)
}

#count genes in the after gene column
count_gene_downstream = function(result,distance_min,distance_max){
  count_up_minus_strand = 0
  count_down_plus_strand = 0
  for(i in 1:length(result$TE_id)){
    downstream_gene_end = result$after_end[i]
    if(is.nan(downstream_gene_end) == FALSE & result$TE_strand[i] == '+' & result$Down_Genestart.TEend[i] > distance_min & result$Down_Genestart.TEend[i] < distance_max){
      count_down_plus_strand = count_down_plus_strand + 1
    }
    if(is.nan(downstream_gene_end) == FALSE & result$TE_strand[i] == '-' & result$Down_Genestart.TEend[i] > distance_min & result$Down_Genestart.TEend[i] < distance_max){
      count_up_minus_strand = count_up_minus_strand + 1
    }
  }
  count = c(count_up_minus_strand,count_down_plus_strand)
  return(count)
}

c1 = count_gene_upstream(result_file,0,500)
c2 = count_gene_upstream(result_file,501,1000)
c3 = count_gene_upstream(result_file,1001,2000)
c4 = count_gene_upstream(result_file,2001,max(result_file$end))

c9 = count_gene_downstream(result_file,0,500)
c10 = count_gene_downstream(result_file,501,1000)
c11 = count_gene_downstream(result_file,1001,2000)
c12 = count_gene_downstream(result_file,2001,max(result_file$end))

pdf(opt$pdf)
#ggplot graph : "Number of LTR with an adjacent upstream gene"
distance = c(rep("0-0.5kbp" , 2), rep("0.5-1kbp" , 2), rep("1-2kbp", 2), rep("2kbp+", 2))
strand = rep(c("+","-"),4)
value_up = c(c1[1],c9[1],c2[1],c10[1],c3[1],c11[1],c4[1],c12[1])
data = data.frame(distance,strand,value_up)
plot1 = ggplot(data, aes(fill=strand, y=value_up, x=distance )) + geom_bar(position="dodge", stat="identity") + ylab("number of LTR") + xlab("distance to the gene") + ggtitle("Number of LTR with an adjacent upstream gene") + theme(plot.title = element_text(hjust = 0.5))
print(plot1)

#ggplot graph : "Number of LTR with an adjacent downstream gene"
distance = c(rep("0-0.5kbp" , 2), rep("0.5-1kbp" , 2), rep("1-2kbp", 2), rep("2kbp+", 2))
strand = rep(c("+","-"),4)
value_down = c(c1[2],c9[2],c2[2],c10[2],c3[2],c11[2],c4[2],c12[2])
data = data.frame(distance,strand,value_down)
plot2 = ggplot(data, aes(fill=strand, y=value_down, x=distance )) + geom_bar(position="dodge", stat="identity") + ggtitle("Mandshurica") + ylab("number of LTR") + xlab("distance to the gene") + theme(plot.title = element_text(hjust = 0.5))
print(plot2)
invisible(dev.off())
print("Plots generated")

#write result of counting in file
df = data.frame(col1 = c(c1[1],c2[1],c3[1],c4[1]), col2 = c(c9[1],c10[1],c11[1],c12[1]), col3 = c(c9[2],c10[2],c11[2],c12[2]), col4 = c(c1[2],c2[2],c3[2],c4[2]))
colnames(df) = c("Upstream overlap(+)","Upstream overlap(-)","Downstream_overlap(+)","Downstream_overlap(+)")
row.names(df) = c("0-500bp","501-1000bp","1001-2000bp","2000+bp")
write.table(df, opt$out, row.names = TRUE, sep="\t", col.names=NA)
print("Table written")

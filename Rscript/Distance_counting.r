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
              help="dataset file name i.e. 'result/output_TE.tsv'", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="result/",
              help="Path to output [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# no option given to script
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file). ", call.=FALSE)
}

#load result file
result_file = read.csv(file = opt$file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#count genes in the before gene column
count_gene_upstream = function(result, distance_min, distance_max){
  count = 0
  for(i in 1:length(result$TE_id)){
    if (result$relationship[i] == "upstream"){
      if (result$TE_strand[i] == "+" & result$Gene_strand[i] == "+"){
        distance = abs(result$Gene_start[i] - result$TE_end[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }else if (result$TE_strand[i] == "-" & result$Gene_strand[i] == "+"){
        distance = abs(result$TE_start[i] - result$Gene_start[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }else if (result$TE_strand[i] == "-" & result$Gene_strand[i] == "-"){
        distance = abs(result$Gene_start[i] - result$TE_end[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }else if (result$TE_strand[i] == "+" & result$Gene_strand[i] == "-"){
        distance = abs(result$Gene_start[i] - result$TE_start[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }
    }
  }
  return(count)
}

#count genes in the after gene column
count_gene_downstream = function(result,distance_min,distance_max){
  count = 0
  for(i in 1:length(result$TE_id)){
    if (result$relationship[i] == "downstream"){
      if (result$TE_strand[i] == "+" & result$Gene_strand[i] == "+"){
        distance = abs(result$Gene_end[i] - result$TE_start[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }else if (result$TE_strand[i] == "-" & result$Gene_strand[i] == "+"){
        distance = abs(result$Gene_end[i] - result$TE_end[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }else if (result$TE_strand[i] == "-" & result$Gene_strand[i] == "-"){
        distance = abs(result$TE_start[i] - result$Gene_end[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }else if (result$TE_strand[i] == "+" & result$Gene_strand[i] == "-"){
        distance = abs(result$TE_end[i] - result$Gene_end[i])
        if ( distance_min <= distance & distance <= distance_max){
          count = count + 1
        }
      }
    }
  }
  return(count)
}

c1 = count_gene_upstream(result_file,0,500)
c2 = count_gene_upstream(result_file,501,1000)
c3 = count_gene_upstream(result_file,1001,2000)
c4 = count_gene_upstream(result_file,2001,max(result_file$Gene_end))

c9 = count_gene_downstream(result_file,0,500)
c10 = count_gene_downstream(result_file,501,1000)
c11 = count_gene_downstream(result_file,1001,2000)
c12 = count_gene_downstream(result_file,2001,max(result_file$Gene_end))

print("Create distance plots.")
#ggplot graph : "Number of TE with an adjacent upstream gene"
distance = c("0-0.5kbp", "0.5-1kbp", "1-2kbp", "2kbp+")
value_up = c(c1, c2, c3, c4)
data = data.frame(distance,value_up)
g1 = ggplot(data, aes(y=value_up, x=distance )) + 
  geom_bar(position="dodge", stat="identity") + 
  ylab("number of TE") + 
  xlab("distance to the gene") + 
  ggtitle("Number of TE with an adjacent upstream gene") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(opt$out, "TEs_upstream.png"), plot = g1, width = 8, height = 8)


#ggplot graph : "Number of TE with an adjacent downstream gene"
distance = c("0-0.5kbp", "0.5-1kbp", "1-2kbp", "2kbp+")
value_down = c(c9, c10, c11, c12)
data = data.frame(distance,value_down)
g2 = ggplot(data, aes(y=value_down, x=distance )) + 
  geom_bar(position="dodge", stat="identity") + 
  ggtitle("Number of TE with an adjacent downstream gene") + 
  ylab("number of TE") + 
  xlab("distance to the gene") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(opt$out, "TEs_downstream.png"), plot = g2, width = 8, height = 8)

#write result of counting in file
df = data.frame(col1 = c(c1,c2,c3,c4), col2 = c(c9,c10,c11,c12))
colnames(df) = c("Upstream", "Downstream")
row.names(df) = c("0-500bp","501-1000bp","1001-2000bp","2000+bp")
write.table(df, paste0(opt$out, "distance_TE_results.csv"), row.names = TRUE, sep=",", col.names=NA)
print("Table written")

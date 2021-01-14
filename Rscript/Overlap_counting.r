if(!require("optparse")){
    install.packages("optparse")
    library(optparse)
}
if(!require("ggplot2")){
    install.packages("ggplot2")
    library(ggplot2)
}
if(!require("stringr")){
    install.packages("stringr")
    library(stringr)
}

# script options
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name i.e. 'result/output_TE.tsv'", metavar="character"),
  make_option(c("-p", "--pdf"), type="character", default="result/overlap_TE_results.pdf",
              help="output filename name (PDF) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="result/overlap_TE_results.csv",
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

#count genes in the upstream overlap column
count_overlap_upstream = function(result){
count_up_plus_strand = 0
count_down_minus_strand = 0
for(i in 1:length(result$TE_id)){
    a = str_remove_all(result$upstream_overlap_feature[i], "[\"',\\[\\]]")
    features = strsplit(a, "\\s+")[[1]]
    TE_strand = result$TE_strand
    if(result$upstream_overlap[i] != '[]'){
        for(j in 1:length(features)){
            if(features[j] == 'exon' & TE_strand[i] == '+'){
                count_up_plus_strand = count_up_plus_strand + 1
                break
                }
            if(features[j] == 'exon' & TE_strand[i] == '-'){
                count_down_minus_strand = count_down_minus_strand + 1
                break
                }
            
            }
        }
    }
    count = c(count_up_plus_strand,count_down_minus_strand)
    return(count)
}

#count genes in the downstream overlap column
count_overlap_downstream = function(result,sens){
count_up_minus_strand = 0
count_down_plus_strand = 0
for(i in 1:length(result$TE_id)){
    a = str_remove_all(result$downstream_overlap_feature[i], "[\"',\\[\\]]")
    features = strsplit(a, "\\s+")[[1]]
    TE_strand = result$TE_strand
    if(result$downstream_overlap[i] != '[]'){
        for(j in 1:length(features)){
            if(features[j] == 'exon' & TE_strand[i] == '+'){
                count_down_plus_strand = count_down_plus_strand + 1
                break
                }
            if(features[j] == 'exon' & TE_strand[i] == '-'){
                count_up_minus_strand = count_up_minus_strand + 1
                break
                }
            
            }
        }
    }
    count = c(count_up_minus_strand,count_down_plus_strand)
    return(count)
}

c1 = count_overlap_upstream(result_file)
c2 = count_overlap_downstream(result_file)

#create graph
sens = c(rep("Upstream Overlap", 2), rep("Downstream Overlap", 2))
strand = rep(c("+","-"),2)
value = c(c1[1],c2[1],c1[2],c2[2])
data = data.frame(sens,value)
pdf(opt$pdf)
plot = ggplot(data, aes(fill=strand, y=value, x=sens )) + geom_bar(position="dodge", stat="identity") + ggtitle("Dataset") + ylab("Number of TE with an overlapping gene") + xlab("") + theme(plot.title = element_text(hjust = 0.5))
print(plot)
invisible(dev.off())
print("INFO - Plot created")

#write result of counting in file
df = data.frame(col1 = c1[1], col2 = c2[1], col3 = c1[2], col4 = c2[2])
colnames(df) = c("Upstream overlap(+)","Upstream overlap(-)","Downstream_overlap(+)","Downstream_overlap(+)")
write.table(df, opt$out,row.names = FALSE, sep="\t")
print("INFO - Table written")

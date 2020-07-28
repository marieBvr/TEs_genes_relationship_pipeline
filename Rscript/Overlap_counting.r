library(ggplot2)
library(stringr)

#load result file
result_file = read.csv(file = 'ResultFile_LTR/ResultFile_LTR_Stella.tsv', sep = '\t', header = TRUE)

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
plot = ggplot(data, aes(fill=strand, y=value, x=sens )) + geom_bar(position="dodge", stat="identity") + ggtitle("Sibirica") + ylab("number of TE with an overlapping gene") + xlab("") + theme(plot.title = element_text(hjust = 0.5))
print(plot)

#write result of counting in file
df = data.frame(col1 = c1[1], col2 = c2[1], col3 = c1[2], col4 = c2[2])
colnames(df) = c("Upstream overlap(+)","Upstream overlap(-)","Downstream_overlap(+)","Downstream_overlap(+)")
write.table(df,"overlap_count_result.csv",row.names = FALSE, sep="\t")

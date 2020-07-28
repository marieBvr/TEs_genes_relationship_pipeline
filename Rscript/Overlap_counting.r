library(ggplot2)
library(stringr)


result_file = read.csv(file = 'ResultFile_TE/ResultFile_Sibirica.tsv', sep = '\t', header = TRUE)
count_overlap_upstream = function(result,sens){
count = 0
for(i in 1:length(result$TE_id)){
    a = str_remove_all(result$upstream_overlap_feature[i], "[\"',\\[\\]]")
    features = strsplit(a, "\\s+")[[1]]
    b = str_remove_all(result$upstream_overlap_strand[i], "[\"',\\[\\]]")
    strands = strsplit(b, "\\s+")[[1]]
    if(result$upstream_overlap[i] != '[]'){
        for(j in 1:length(features)){
            if(features[j] == 'exon' & strands[j] == sens){
                count = count + 1
                break
                }
            
            }
        }
    }
    print(count)
    return(count)
}

count_overlap_downstream = function(result,sens){
count = 0
for(i in 1:length(result$TE_id)){
    a = str_remove_all(result$downstream_overlap_feature[i], "[\"',\\[\\]]")
    features = strsplit(a, "\\s+")[[1]]
    b = str_remove_all(result$downstream_overlap_strand[i], "[\"',\\[\\]]")
    strands = strsplit(b, "\\s+")[[1]]
    if(result$downstream_overlap[i] != '[]'){
        for(j in 1:length(features)){
            if(features[j] == 'exon' & strands[j] == sens){
                count = count + 1
                break
                }
            
            }
        }
    }
    print(count)
    return(count)
}

c1 = count_overlap_upstream(result_file,'+')
c2 = count_overlap_upstream(result_file,'-')

c3 = count_overlap_downstream(result_file,'+')
c4 = count_overlap_downstream(result_file,'-')

count = c(c1,c2,c3,c4)

sens = c(rep("Upstream Overlap", 2), rep("Downstream Overlap", 2))
strand = rep(c("+","-"),2)
value = count
data = data.frame(sens,value)
plot = ggplot(data, aes(fill=strand, y=value, x=sens )) + geom_bar(position="dodge", stat="identity") + ggtitle("Sibirica") + ylab("number of TE with an overlapping gene") + xlab("") + theme(plot.title = element_text(hjust = 0.5))
print(plot)


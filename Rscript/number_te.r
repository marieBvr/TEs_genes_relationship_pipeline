library(ggplot2)
#load result file
data = read.table(file = 'ResultFile_Msa.tsv', sep = '\t', header = TRUE)
df = data[with(data, order(data$TE_Type)), ]


find_element = function(df){
  elements = c()
  element = df$TE_Type[1] #print the first row of the TE_Type column
  print(element)
  elements <- c(elements,element)
  
  for(i in 1:length(df$TE_Type)){
    if(df$TE_Type[i]!=element){
      element = df$TE_Type[i]
      
      elements <- c(elements,element)
    }
  } 
  print((elements))
  return (elements)
}


number_of_element = function(df, element){
  compter = 0
  for(i in 1:length(df$TE_Type)){
      if(df$TE_Type[i]==element){
        compter = compter + 1
      }
    }
    return(compter)
}


#main 
 list_of_element = find_element(df)
 number = c()
 for( i in list_of_element){
   a = number_of_element(df, i)
   number <- c(number,a)
 }
 print(number)


 #graph 
 dat <- data.frame(x=list_of_element, y=number)
 barplot(dat$y, names.arg=dat$x,
         main = "Number of transposon for each type of TE ",
         ylab="number", 
         xlab="type",
         col=list_of_element)








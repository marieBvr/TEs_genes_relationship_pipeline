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
  make_option(c("-o", "--out"), type="character", default="Resting_result/", 
              help="output directory name [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# no option given to script
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file). ", call.=FALSE)
}

#load result file
filename = opt$file
data = read.table(file = filename, sep = '\t', header = TRUE)
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
  #print(elements)
  return(elements)
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
#print(number)


#graph 
dat <- data.frame(x=list_of_element, y=number)
# Open a pdf file
pdf("count_TE_transposons.pdf")
barplot(dat$y, names.arg=dat$x,
        main = "Number of transposon for each type of TE ",
        ylab="number", 
        xlab="type",
        col=list_of_element)
dev.off()
print("Number of elements calculated, plot generated.")

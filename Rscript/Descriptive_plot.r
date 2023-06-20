if(!require("optparse")){
    install.packages("optparse")
    library(optparse)
}
if(!require("ggplot2")){
    install.packages("ggplot2")
    library(ggplot2)
}

# script options
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name i.e. 'result/output_TE.tsv'", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="result/",
              help="output path [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# no option given to script
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file). ", call.=FALSE)
}

#load result file
file = read.csv(file = opt$file, sep = '\t', header = TRUE)

# Create plots
print("Create 'Relationships' plot.")
g1 <- ggplot(data = file)+
  geom_bar(aes(x = relationship, fill = relationship))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = "Relationships")+
  guides(fill="none")

ggsave(paste0(opt$out, "boxplot_relationships.png"), plot = g1, width = 8, height = 8)

print("Create 'Amount of TEs by gene' plot.")
g2 <- ggplot(data = file, aes(x = Gene_id, fill = Type)) +
    geom_bar(position = "fill") + ylab("proportion") +
    stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=0.5), colour="white")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(title = "Amount of TEs by gene")

ggsave(paste0(opt$out, "boxplot_TE_by_Gene.png"), plot = g2, width = 8, height = 8)

print("Create 'Amount of TEs by Chr' plot.")
g3 <- ggplot(data = file, aes(x = Chr, fill = Type)) +
    geom_bar(position = "fill") + ylab("proportion") +
    stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=0.5), colour="white")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(title = "Amount of TEs by Chr")

ggsave(paste0(opt$out, "boxplot_TE_by_Chr.png"), plot = g3, width = 8, height = 8)

print("Create 'Relationships type of TEs' plot.")
g4 <- ggplot(data = file, aes(x = Type, fill = relationship)) +
    geom_bar(position = "fill") + ylab("proportion") +
    stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=0.5), colour="white")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(title = "Relationships type of TEs")

ggsave(paste0(opt$out, "boxplot_relationships_by_TEs.png"), plot = g4, width = 8, height = 8)

print("Create 'Frequency of TEs' plot.")
g5 <- ggplot(data = file)+
  geom_bar(aes(x = Type, fill = Type))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = "Frequency of each TEs")+
  guides(fill="none")

ggsave(paste0(opt$out, "boxplot_frequency_TEs.png"), plot = g5, width = 8, height = 8)

print("Done.")
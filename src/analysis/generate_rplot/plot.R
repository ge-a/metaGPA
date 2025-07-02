library(ggplot2)
library(gggenes)
library(RColorBrewer)
library(dbplyr)
library(stringr)
#library(tidyverse)

#setwd('/Users/ettwiller/Misc/projects/metaGPA/dZ_DNA/contigs')
#file_name <- "test_PF05272.14.png"

args <- commandArgs(trailingOnly=TRUE)
filein <- args[1]
file_name <- args[2]

mutation <- read.table(filein, header=FALSE, sep="")
colnames(mutation) <- c("molecule","gene","start","end","strand", "orientation","selection", "color");
#typ = unique(mutation$type)

mutation$ratio <- str_extract(mutation$molecule, "ratio_[.,0-9]+")
mutation$ratio <- as.numeric(str_extract(mutation$ratio, '\\d+([.,]\\d+)?'))

pal <- c(
  "PF03602.18" = "red",
  "PF00145.20" = "orange", 
  "PF16631.8" = "yellow", 
  "PF17728.4" = "forestgreen",
  "PF09517.13" ="NA"
)
test <- data.frame(pal)

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n)

contigs <- length(unique(mutation$molecule)) 
geneheight_ratio = contigs/3
geneheight <- (contigs/geneheight_ratio)
imageheight <- (contigs/5)
arrowh <- ((geneheight/9))

#mutation$selection_reorder = factor(mutation$selection, levels=c('DEPLETED','SELECTED'))
mutation$selection_reorder = factor(mutation$selection, levels=c('SELECTED','DEPLETED'))

t<- ggplot(mutation, aes(xmin = start, xmax = end, y = reorder(molecule,ratio), fill = color, forward = orientation)) +
  geom_gene_arrow(arrowhead_height = unit(geneheight, "mm"), arrowhead_width = unit(arrowh, "mm"),arrow_body_height = unit(geneheight, "mm")) +
  facet_grid(selection_reorder~., scales = "free", space='free') +
  #facet_wrap(selection~ .) +
  #scale_fill_manual(values = col) +
  #scale_fill_manual(values = pal,limits = names(pal),na.value=NA) 
  scale_fill_identity()+
  theme_genes()+
  geom_gene_label(label=mutation$gene, size =5)+
  theme(
        #panel.background = element_rect(fill = 'white', color = 'grey40', linewidth =0.2),
        panel.grid.minor = element_line(color = 'grey80',size =0.2),
        panel.grid.major = element_line(color = 'grey60', size =0.2),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),)

ggsave(file_name, t, width=15, height=imageheight, limitsize = FALSE)
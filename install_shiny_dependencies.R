#install dependencies
#run each code chunk separately 

##########
#chunck 1#
##########
for (p in c("shiny", "DT", "htmltools", "ggplot2","data.table","tidyverse","httr","jsonlite","curl","RCurl","ggbeeswarm"))
{
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, dependencies = TRUE)
}

##########
#chunck 2#
##########
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(BiocManager)

##########
#chunck 3#
##########
BiocManager::install("DESeq2")

##########
#chunck 4#
##########
BiocManager::install("biomaRt")

##########
#chunck 5#
##########
BiocManager::install("EnhancedVolcano")
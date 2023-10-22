


library(vegan)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(scales)
library(stringi)
library(MASS)
#library(stargazer)
library(reporttools)
library(epitools)
library(gdata)
library(car)
library(plyr)
library(dplyr)
library(data.table)
library(tibble)
library(psych)
library(tidyr)
library(janitor)
library(psych)
library(plotrix)
library(slopegraph)
library(Lock5Data)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(treemap)
library (treemapify)
library(ggraph)
library(igraph)



###The code, uses the following "INPUT" files:

# CARD_read_count_specific.tsv
# ResPipe_CARD-3.0.3.meta.tsv
# CARD_lateral_coverage_specific.tsv
# CARD_read_lengths_specific.tsv
# Metagenomics_Metadata.csv
# Infection_Data_For_Bayesian_Model.csv
# bracken_combined_reads.tsv

# to produce the corrected resistance gene counts, plus a matrix that links each resistance gene and antibiotic based on the "Confers_Resistance_to_Antibiotic" relationship ontology term.
    # in the matrix, 1 == the gene is associated with clear experimental evidence of elevated MIC for that antibiotic but the "Confers_Resistance_to_Antibiotic" relationship ontology term is missing
    # in the matrix, 2 == the gene is associated with demonstrably elevated MIC for that antibiotic and is known to confer or contribute to clinically relevant resistance to that antibiotic ("Confers_Resistance_to_Antibiotic" relationship ontology term is present)

#OUTPUT FILES:

# Corrected_Counts.csv
# AB_Matrix_1_or_2.csv

# AMR_DEF.csv
# AMR_ALL.csv

# And to produce the final dataset for the bayesian modelling (i.e. "OUTPUT" file)

# Dataset_For_Bayesian_Model.csv

# The code for the non-metric multidimensional scaling (NMDS) ordination method is also presented (see supplementary results for the validation of the pooling using 30-sample pools)


##################PRODUCE THE DATASET WITH CORRECTED GENE COUNTS#############################

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix <<- prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }


#CARD_read_count_specific<-read.csv("inputs/CARD_read_count_specific.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
#Gene_Lenghts<-read.csv("inputs/ResPipe_CARD-3.0.3.meta.tsv", header=T,check.names = F, row.names = 1, sep = "\t")
#Spec_Lat_Cov<-read.csv("inputs/CARD_lateral_coverage_specific.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
#Read_Lenghts<-read.csv("inputs/CARD_read_lengths_specific.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
#metadata = read.csv(file = "inputs/Metagenomics_Metadata.csv", header = T, row.names = 1, check.names = F)

CARD_read_count_specific<<-read.csv(paste(pfix, parameters["readcounts", 2], sep="/"), sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
Gene_Lenghts<<-read.csv(paste(pfix, parameters["genelengths", 2], sep="/"), header=T,check.names = F, row.names = 1, sep = "\t")
Spec_Lat_Cov<<-read.csv(paste(pfix, parameters["lateralcoverage", 2], sep="/"), sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
Read_Lenghts<<-read.csv(paste(pfix, parameters["readlengths", 2], sep="/"), sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
metadata <<- read.csv(file = paste(pfix, parameters["metadata", 2], sep="/"), header = T, row.names = 1, check.names = F)
RDSdir <<- paste(pfix, parameters["RDS", 2], sep="/")
}


run <- function() {}

output <- function(outputfile) {
MT<-Gene_Lenghts
GL<-subset(MT, select=("SeqLength"))

SLC<-Spec_Lat_Cov

ARL<-Read_Lenghts



memory.limit(size = 1000000000)
names(SLC) = gsub(pattern = "CARD_", replacement = "", x = names(SLC))
names(ARL) = gsub(pattern = "CARD_", replacement = "", x = names(ARL))
names(CARD_read_count_specific) = gsub(pattern = "CARD_", replacement = "", x = names(CARD_read_count_specific))

identical(rownames(GL), rownames(ARL))
identical(rownames(ARL), rownames(SLC))
identical(rownames(SLC), rownames(CARD_read_count_specific))
identical(colnames(SLC), colnames(ARL))
identical(colnames(SLC), colnames(CARD_read_count_specific))

GSL<-mapply("*", SLC,GL)
GLARL=GSL/ARL 
GLARL[is.na(GLARL)] <- 0
Corr_Count=CARD_read_count_specific/GLARL 
Corr_Count[is.na(Corr_Count)] <- 0


##########################MAP CORRECTED COUNTS TO AROs################################################################

newnamesd1<-names(Corr_Count)

#attach(metadata)

#fix(metadata)
z = metadata[match(newnamesd1, rownames(metadata)),]


#REMOVE 24 FOLLOW-UP SAMPLES FROM PARALLEL STUDY


#data(Corr_Count)
dropList <- c("NCS-ST-0047",
              "NCS-ST-0049",
              "NCS-ST-0051",
              "NCS-ST-0065",
              "NCS-ST-0073",
              "NCS-ST-0079",
              "NCS-ST-0110",
              "NCS-ST-0175",
              "NCS-ST-0176",
              "NCS-ST-0191",
              "NCS-ST-0247",
              "NCS-ST-0251",
              "NCS-ST-0257",
              "NCS-ST-0380",
              "NCS-ST-0389",
              "NCS-ST-0413",
              "NCS-ST-0416",
              "NCS-ST-0417",
              "NCS-ST-0423",
              "NCS-ST-0426",
              "NCS-ST-0429",
              "NCS-ST-0438",
              "NCS-ST-0454",
              "NCS-ST-0457")

Int_1 <- Corr_Count[, !colnames(Corr_Count) %in% dropList]
Corr_Count<-Int_1

#####merge gene data with ARO labels#######################


df<-MT
df[,2:12]<- list(NULL)
labels_RG<-df

All_Data <- merge(Corr_Count, labels_RG, by=0, all=FALSE)


#remove Row names

All_Data$Row.names<-NULL

## move ARO accession number column to first position

Int_1 <- All_Data %>%
  select(ARO_accession, everything()) 

####a few aro accession numbers map up to more than one R_XXX UID. Hence we need a dataset that
#### sums counts across all rows that map to the same aro number


Int_2<-aggregate(.~ARO_accession, data=Int_1, FUN=sum) 

##convert aroaccession to row names

row.names(Int_2) <- Int_2[,1]
Int_2["ARO_accession"]<-NULL 


##remove target genes where the sum across ALL SAMPLES AND POOLS is zero. 


total_col<-apply(Int_2[,], 1, sum)
Int_3<-as.data.frame(cbind(Int_2,total_col))
Int_4<-Int_3[Int_3$total_col!=0, ]

nrow(Int_4) ##Total number of ARO accession IDs detected in at least one of the samples

ARO_Counts<-Int_4
ARO_Counts$total_col<-NULL

#ARCHIVE ARO DATA

write.csv(ARO_Counts, outputfile) 
saveRDS(ARO_Counts, paste(RDSdir, "ARO_Counts.rds", sep="/"))

}


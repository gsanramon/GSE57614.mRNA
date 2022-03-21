### I. Get annotations of proteins in Endoplasmic Reticulum (ER), specifically located in the membrane. 

#Load ER list
ER <- read.delim('ER.txt')
dim(ER) #10324 annotations and 14 columns
head(ER)
#Find unique receptor proteins by searching for the phrase 'endoplasmic reticulum membrane' under column 'GO.NAME'
ER.membrane <- ER[which(ER$GO.NAME == 'endoplasmic reticulum membrane'),]
length(unique(ER.membrane$SYMBOL)) #This is the number of proteins in the ER membrane
head(ER.membrane)


### II. Find genes that are highly expressed or repressed in M1 macrophages using RNA data that is publicly available from NCBI GEO website.

## M0 vs M1
#Load the file for M0 vs M1 (24hrs)
M0vsM1.24hr <- read.delim('GSE57614.M0vsM1.24hr.txt')
head(M0vsM1.24hr)
#Extract < 5% FDR
M0vsM1.24hr.fdr5 <- M0vsM1.24hr[which(M0vsM1.24hr$adj.P.Val < 0.05),]
#dim(M0vsM1.24hr.fdr5) 

#Get the significant genes that are unique minus the blanks ''
M0M1 <- unique(M0vsM1.24hr.fdr5$Gene.symbol[-which(M0vsM1.24hr.fdr5$Gene.symbol=='')])
length(M0M1) # 579 genes are significant

##M2 vs M1
#Load the file for M2 vs M1 (24hrs)
M2vsM1.24hr <- read.delim('GSE57614.M2vsM1.24hr.txt')
head(M2vsM1.24hr)
#Extract < 5% FDR
M2vsM1.24hr.fdr5 <- M2vsM1.24hr[which(M2vsM1.24hr$adj.P.Val < 0.05),]
#dim(M2vsM1.24hr.fdr5) 

#Get the significant genes that are unique minus the blanks ''
M2M1 <- unique(M2vsM1.24hr.fdr5$Gene.symbol[-which(M2vsM1.24hr.fdr5$Gene.symbol=='')])
length(M2M1) # 926 genes are significant

#Find the common genes in the two lists:
common2 <- Reduce(intersect, list(M2M1, M0M1))
length(common2) #There are 429 genes common to M2vsM1 and M0vsM1
#Remove hash to create a file
#write.csv(common2, 'common2.csv')

#show all 429 genes common to M2vsM1 and M0vsM1 lists
#common2

#Find intersection between variables 'common2' and 'ER.membrane'
ER.membrane.X.common2 <- Reduce(intersect, list(ER.membrane$SYMBOL, common2))
length(ER.membrane.X.common2) #Number of DEGs in ER membranev

#List of DEGs in ER membrane 
ER.membrane.X.common2


#Get the data for the 25 ER membrane genes from **M2vsM1** RNA file
data.M2vsM1.DEG.ERmem <- dplyr::filter(M2vsM1.24hr.fdr5, M2vsM1.24hr.fdr5$Gene.symbol %in% ER.membrane.X.common2)
#sort according to logFC value
sort.M2vsM1.DEG.ERmem <- data.M2vsM1.DEG.ERmem[order(data.M2vsM1.DEG.ERmem$logFC),]
head(sort.M2vsM1.DEG.ERmem)
#To create a file, remove hash below
#write.csv(sort.M2vsM1.DEG.ERmem, 'sort.M2vsM1.DEG.ERmem.csv')

#Get the data for the 25 ER membrane genes from **M0vsM1** RNA file
data.M0vsM1.DEG.ERmem <- dplyr::filter(M0vsM1.24hr.fdr5, M0vsM1.24hr.fdr5$Gene.symbol %in% ER.membrane.X.common2)
#sort according to logFC value
sort.M0vsM1.DEG.ERmem <- data.M0vsM1.DEG.ERmem[order(data.M0vsM1.DEG.ERmem$logFC),]
head(sort.M0vsM1.DEG.ERmem)
#To create a file, remove hash below
#write.csv(sort.M0vsM1.DEG.ERmem, 'sort.M0vsM1.DEG.ERmem.csv')

##Load both antibody files and cross reference the 25 DEGs.**
  
#  1. Using Ab_Chris_refined.csv file:

#load Ab data
Ab <- read.csv('Ab_Chris_refined.csv')
dim(Ab)

#The gene/s below is differentially expressed and located in the ER membrane.

#cross reference the 25 DEGs to the antibody list
ERmem.com.Ab <- Reduce(intersect, list(ER.membrane.X.common2, Ab$network_name))
#show available antibody
ERmem.com.Ab

#2. Using phospho_Chris_refined.csv file:

#load data for phospho Ab (this table doesn't have header)
pAb <- read.csv('phospho_Chris_refined.csv', header=FALSE)
dim(pAb)
#add column names
colnames(pAb) <- c('phospho.name', 'gene', 'phosphorylation')

#The gene/s below is differentially expressed and located in the ER membrane.

#cross reference the 25 DEGs to the phospho antibody list
ERmem.com.pAb <- Reduce(intersect, list(ER.membrane.X.common2, pAb$gene))
#show available antibody
ERmem.com.pAb





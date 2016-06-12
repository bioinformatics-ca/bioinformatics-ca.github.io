######## integrated assignment ######
####
#install required R and bioconductor packages
source("http://www.Bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

biocLite(c("Biobase","limma", "pheatmap",  "RColorBrewer"))

#load Required libraries
library("Biobase")
library("limma")
library("pheatmap")
library("RColorBrewer")
library("affy")
library("biomaRt")

## set working directory
setwd("~/Dropbox (Bader Lab)/Veronique Voisin's files/CBW/new_integrated_assignment");


### read data
matrix <- read.delim("GSE39491_series_matrix.txt",  quote=NULL, stringsAsFactors=FALSE)
head(matrix)

# id column have some quotes that I need to remove and the sample ids have some weird X and. that I need to remove as well.
colnames(matrix) <- gsub("X\\.", "", colnames(matrix))
colnames(matrix) <- gsub("\\.", "",  colnames(matrix))
rownames(matrix) <- gsub("\"", "", matrix$ID_REF)
matrix$ID_REF <- gsub("\"", "", matrix$ID_REF)

metadata <- read.csv("metadata.csv",  stringsAsFactors=FALSE)
head(metadata)

length(which(metadata$sample %in% colnames(matrix))) #120

NE_group <- metadata$sample[which(metadata$type == "NE")]
BE_group <- metadata$sample[which(metadata$type == "BE")]

#################################################
# BE vs NE
#################################################

#create the data matrix for limma with the expression value
minimalSet <- ExpressionSet(assayData=log2(as.matrix(matrix[ , -c(1)])))

classes <- factor(metadata[,"type"])

modelDesign <- model.matrix(~ 0 + classes)

colnames(modelDesign)

fit <- lmFit(minimalSet, modelDesign) 

#(2) contrast fitting
contrastnm <- c("classesBE-classesNE") #contrast between a pair of classes
contrast.matrix <- makeContrasts(contrasts=contrastnm, levels=modelDesign)

fit1 <- contrasts.fit(fit, contrast.matrix)

#(3) eBayes fitting
fit2 <- eBayes(fit1)

#(4.1) topfit based on F-statistic
#H1: At least one of the contrasts is significantly differetial.  
topfit <- topTable(fit2, number=nrow(minimalSet), adjust="BH")
head(topfit)
## get the gene id from biomart
## get annotation from biomart (2)
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
listAttributes(mart)[grep("*affy*", listAttributes(mart))]
genes = getBM(attributes = c('affy_hg_u133a', 'hgnc_symbol', 'description'), filters='affy_hg_u133a', values=row.names(topfit), mart=mart);
genes[genes==""] = NA;
head(genes)
names(genes)

sum(duplicated(genes$hgnc_symbol))
sum(duplicated(genes$description))
genes$description = gsub("\\[Source.*", "", genes$description);
sum(is.na(genes$hgnc_symbol))

### reduce at the gene level
##duplicated probes
##- one probe ID might correspond to more than one gene that could be in different chromosomes!
rm(Gene)

## annnotated normalized data and remove duplicated gene entries
m = match(row.names(topfit), genes$affy_hg_u133a)

##- collapse gene names for the duplicated probe ID
dupProbe <- unique( genes$affy_hg_u133a[duplicated(genes$affy_hg_u133a)] )
Gene <-  genes$hgnc_symbol
for (i in 1:length(dupProbe)) {
   indx <- (genes$affy_hg_u133a == dupProbe[i])
   Gene[indx] <- paste0(genes$hgnc_symbol[indx], collapse = "/")
   }
genes$upd_symbol <- Gene

topfit$symbol <- genes$upd_symbol [m]
topfit$description <- genes$description [m]
head(topfit)

### nerge topfit with the normalized data
topfit$ID<- rownames(topfit)
topfit_expression <- merge(topfit, matrix, by.x='ID', by.y='ID_REF')
head(topfit_expression, n=2)

topfit_expression <- topfit_expression[order(abs(topfit_expression$t), decreasing=TRUE)  , ]

topfit_expression <- topfit_expression[ !is.na(topfit_expression$symbol) & !duplicated(topfit_expression$symbol) , ]

topfit_expression <- topfit_expression[order(topfit_expression$t, decreasing=TRUE)  , ]

write.csv(topfit_expression, "topfit_expression.csv")
#create the output matrix
ranks <- topfit_expression[, c("symbol", "t")]

write.table(ranks,"BEvsNE_ranks.rnk",col.name=FALSE,sep="\t",row.names=FALSE,quote=FALSE)

grep(NE_group, colnames(topfit_expression)[NE_group] )
#output optional expresion file to be used by Enrichment Map
EM_expressionFile <- topfit_expression[ , c("symbol", "description",  NE_group, BE_group)] 
colnames(EM_expressionFile)[1] <- "Name"
colnames(EM_expressionFile)[2] <- "Description"
write.table(EM_expressionFile,"BE_vs_NE_expression.txt",col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)


#########################################
# gprofiler files - top genes only; ordered usiong FDR so we can use the gprofiler ordered option


# BE only from BE vs NE
BEvsNE_BE <- topfit_expression[which(topfit_expression$adj.P.Val < 0.00001 & topfit_expression$logFC > 0.585),]
BEvsNE_BE <- BEvsNE_BE[order(BEvsNE_BE$adj.P.Val), ] 

write.table(BEvsNE_BE$symbol,"BEonly_genelist.txt",col.name=FALSE,sep="\t",row.names=FALSE,quote=FALSE)
 
 # NC only from BE vs NE
BEvsNE_NE <- topfit_expression[which(topfit_expression$adj.P.Val < 0.00001 & topfit_expression$logFC < -0.585),]
BEvsNE_NE <- BEvsNE_NE[order(BEvsNE_BE$adj.P.Val), ] 

write.table(BEvsNE_NE$symbol,"NEonly_genelist.txt",col.name=FALSE,sep="\t",row.names=FALSE,quote=FALSE)
 
 
############################################
###################################################################################################################################################### 
## install packages, pull in data
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq")
#BiocManager::install("DESeq2")
# BiocManager::install("limma")
#install.packages('tidyverse')

##load in the libraries
library(DESeq2)
library(tidyverse)
library(limma)


workdir <- "/scratch/Users/tajo5912/peptide/RNA-seq/deseq2/"


##load the bed file
countdata <- read.delim("/scratch/Users/tajo5912/peptide/RNA-seq/counts/featCounts_gene_names_included.txt", 
                        sep="\t", header=TRUE, check.names=FALSE)

#Check data
head(countdata)

# (Re)name Columns
names(countdata) <- c("id", "length", "DP1", "DP2", "DW1", "DW2", 
                      "NP1", "NP2", "NW1", "NW2", 
                      "H2O_DMSO_1", "H2O_DMSO_2", "H2O_DMSO_3", 
                      "H2O_Nutlin_1", "H2O_Nutlin_2", "H2O_Nutlin_3",
                      "Peg_6_DMSO_1", "Peg_6_DMSO_2", "Peg6_DMSO_3",
                      "Peg6_Nutlin_1", "Peg6_Nutlin_2", "Peg6_Nutlin_3")


##Move gene names to new column to remove duplicates

countdata <- mutate(countdata,id=as.character(id))
countdata <- mutate(countdata,gene=sapply(strsplit(countdata$id, split='_', fixed=TRUE),function(x) (x[3])))
countdata <- mutate(countdata,a_type=sapply(strsplit(countdata$id, split='_', fixed=TRUE),function(x) (x[1])))
countdata <- mutate(countdata,a_num=sapply(strsplit(countdata$id, split='_', fixed=TRUE),function(x) (x[2])))
countdata$accession <- paste(countdata$a_type,countdata$a_num, sep='_')
head(countdata)

# remove duplicate genes, keep one with greatest counts/length
# These column numbers will have to be adjusted based on the number of samples you have

countdata$sums <- rowSums( countdata[,3:22] )
countdata$normalized_sums <- (countdata$sums / countdata$length)
head(countdata)

countdata = countdata[order(countdata['gene'],-countdata['normalized_sums']),]
countdata = countdata[!duplicated(countdata[,'gene']),]
countdata = countdata[!duplicated(countdata[,'id']),]
countata <- subset(countdata, select = -c("a_type","a_num"))
head(countdata)

#here I am saving the filtered isoforms
write.table(countdata, paste0(workdir, file="isoform_filtered_counts.txt"), append = FALSE, sep = "\t", row.names=FALSE)

#countdata <- read.csv(paste0(workdir, file="isoform_filtered_counts.csv"))

## Check countdata
head(countdata)

## Remove sums and gene names columns
countdata <- countdata[-c(23:28)]

# Convert column 1 to rownames & remove columns 1 (gene name + accession) and 2 (length)

rownames(countdata) <- countdata[,1]
countdata <- countdata[-c(1,2)]

head(countdata)
dim(countdata)
countdata <- na.omit(countdata)

#Here I am filtering apart my two exper. to run deseq on them separately
#This can be good because if data sets are too difference all the detected variance comes from between datasets rather than between samples of interest (since deseq assumes most things aren't changing)
cols_expt1 <- c(1:8)
cols_expt2 <- c(9:20)
cols_expt2


countdata <- countdata[,cols_expt2]
head(countdata)

# Convert to matrix for DESeq2 analysis
countdata <- as.data.frame(countdata)
countdata <- as.matrix(countdata)



###Limma batch correction
batch <- (c(1,2,3,1,2,3,1,2,3,1,2,3))
#design <- model.matrix(~batch+seq_batch)

countdata_corrected <- removeBatchEffect(countdata, batch=batch)
head(countdata)
head(countdata_corrected)
countdata_corrected <- as.data.frame(countdata_corrected)
countdata_corrected <- round(countdata_corrected,digits=0)
countdata_corrected[countdata_corrected < 0] <- NA
countdata_corrected <- na.omit(countdata_corrected)
head(countdata_corrected)

dim(countdata_corrected)
head(countdata_corrected)

write.table(countdata_corrected, paste0(workdir, file="exp2_limma_corrected_counts.txt"), append = FALSE, sep = "\t")

countdata_corrected <- as.matrix(countdata_corrected)


#conditions table
#expt2
condition <- factor(c(rep("H2O_DMSO", 2), rep("H2O_Nutlin", 2),
                      rep("Peg6_DMSO", 2), rep('Peg6_Nutlin', 2)))
#expt1
#condition <- factor(c(rep("DP", 2), rep("DW", 2), rep("NP", 2), rep("NW", 2)))

condition


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

##############################################################################################################################################
### DESEQ QC
#QC -- Check sizeFactors, should be close to 1 (ideally) -- largely affected by depth (read counts)
sink(paste0(workdir, "size_factors.txt"))
sizeFactors(dds)
sink()

# Plot dispersions
png(paste0(workdir, "qc-dispersions.png", step=""), 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))

png(paste0(workdir, 'hist_reg_log_transform.png'))
hist(assay(rld), col="skyblue", border="slateblue", main="Histogram")
dev.off()

#Examine the effect of log transformation (QC check)
replicate <- 'Nutlin'
place_in_frame <- c(4:5)

png(paste0(workdir, 'log_transform_qc_', replicate, '.png'))
plot(log2( 1 + counts(dds, normalized=TRUE)[ , place_in_frame] ),
     col=rgb(0,0,0,.2), pch=16, cex=0.3 )
dev.off()

# Sample distance heatmap
sampleDists <- dist(t(assay(rld) ) )
png(paste0(workdir,"qc-heatmap-samples.png"), w=1000, h=1000, pointsize=20)
heatmap(as.matrix(sampleDists), key=F, trace="none", main="Sample Distance Matrix")
dev.off()

# Principal components analysis
png(paste0(workdir, "qc-pca.png"), 1000, 1000, pointsize=20)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()


################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################

condx <- "Peg6_Nutlin"   #treatment
condy <- "H2O_Nutlin"     #control


outdir <- paste(workdir, condy, '_vs_', condx, '/', sep='')
dir.create(outdir)

res <- results(dds, contrast = c("condition", condx, condy))
table(res$padj<0.05) 
res <- res[order(res$padj), ]

## Merge with normalized count data
# Notice here as long as we set 'accession_genename' naming format, we'll split that here and save it as two columns
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Accession"
head(resdata)
resdata <- mutate(resdata,Accession=as.character(Accession))
resdata <- mutate(resdata,Gene=sapply(strsplit(resdata$Accession, split='_', fixed=TRUE),function(x) (x[3])))
resdata$Accession <- sub("_[^_]+$", "", resdata$Accession)
resdata <- resdata %>%
  select(Gene, everything())
head(resdata)

### Write results
write.table(resdata, paste0(outdir, file="diffexpr-results_DESeq2.txt"), append = FALSE, sep = "\t", row.names=FALSE)

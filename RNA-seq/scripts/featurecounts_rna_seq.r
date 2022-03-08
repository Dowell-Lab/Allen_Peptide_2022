# R script for featureCounts
# Set up the environment
library(Rsubread)

bamdir <- '/scratch/Shares/dowell/tjones/peptide/mapped/bams/'

input <- c(paste0(bamdir, "DP1_S5_R1_001.sorted.bam"),
           paste0(bamdir, "DP2_S6_R1_001.sorted.bam"),
           paste0(bamdir, "DW1_S1_R1_001.sorted.bam"),
           paste0(bamdir, "DW2_S2_R1_001.sorted.bam"),
           paste0(bamdir, "NP1_S7_R1_001.sorted.bam"),
           paste0(bamdir, "NP2_S8_R1_001.sorted.bam"),
           paste0(bamdir, "NW1_S3_R1_001.sorted.bam"),
           paste0(bamdir, "NW2_S4_R1_001.sorted.bam"),
           paste0(bamdir, "H2O_DMSO_1_S1_R1_001.sorted.bam"),
           paste0(bamdir, "H2O_DMSO_2_S9_R1_001.sorted.bam"),
           paste0(bamdir, "H2O_DMSO_S1_R1_001.sorted.bam"),
           paste0(bamdir, "H2O_Nutlin_1_S2_R1_001.sorted.bam"),
           paste0(bamdir, "H2O_Nutlin_2_S10_R1_001.sorted.bam"),
           paste0(bamdir, "H2O_Nutlin_S2_R1_001.sorted.bam"),
           paste0(bamdir, "Peg_6_DMSO_1_S5_R1_001.sorted.bam"),
           paste0(bamdir, "Peg_6_DMSO_2_S13_R1_001.sorted.bam"),
           paste0(bamdir, "Peg6_DMSO_S7_R1_001.sorted.bam"),
           paste0(bamdir, "Peg6_Nutlin_1_S6_R1_001.sorted.bam"),
           paste0(bamdir, "Peg6_Nutlin_2_S14_R1_001.sorted.bam"),
           paste0(bamdir, "Peg6_Nutlin_S8_R1_001.sorted.bam")
           )

outdir <- "/scratch/Shares/dowell/tjones/peptide/counts/"

# We use a GTF -- this is necessary for featureCounts and more useful for RNA-seq data as it is a little cleaner/faster for counting over exons
hg38gtf <- "/scratch/Shares/dowell/genomes/hg38/hg38_refseq.gtf"

# Read counting using featureCounts
# Sink will save out stdout -- check all of these settings (e.g. isPairedEnd)!
sink(paste0(outdir, "featureCounts_gene_rnaseq.txt"))
fc <- featureCounts(files=input,
                    annot.ext=hg38gtf,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=TRUE,
                    isPairedEnd=FALSE,
                    nthreads=8)
sink()

### Write results -- we'll keep our accession number (GeneID) and length for filtering genes in DESeq2
# Using this output, we can also add gene names (e.g. Name2) using my python script, add_full_genenames.py
write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],
                         fc$counts,stringsAsFactors=FALSE),
            paste0(outdir, file="featureCounts_gene_rnaseq.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)

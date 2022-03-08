library(tidyverse)

workdir <- '/Users/taylorjones/Desktop/Peptide_Project/GSEA_input/'
setwd(workdir)

deseqout <- '/Users/taylorjones/Desktop/Peptide_Project/deseq_expt2_limma_correction/DWvDP/diffexpr-results_DESeq2.txt'
res <- read.csv(deseqout)
head(res)

rnkdf <- tibble(gene = res$Gene,
                rnk = -log(res$pvalue) * sign(res$log2FoldChange)) %>%
  arrange(desc(rnk)) %>% drop_na()
head(rnkdf)

## Write out the table without any additional information
write.table(rnkdf, file = "DW_DP2.rnk", append = FALSE, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")

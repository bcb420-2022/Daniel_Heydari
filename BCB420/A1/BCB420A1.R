if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

if (!requireNamespace("GEOmetadb", quietly = TRUE))
  BiocManager::install("GEOmetadb")
library(GEOmetadb)

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")
library('GEOquery')

if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")
library(edgeR)

# Load in the data
GSE158890 <- getGEO("GSE158890", GSEMatrix =FALSE)

# Print GEO description
head(Meta(GSE158890))

# Print the Platform Title
current_gpl <- names(GPLList(GSE158890))[1]
current_gpl_info <- Meta(getGEO(current_gpl))
current_gpl_info$title

# Print number of GEO datasets, and samples using this technology, respectively
length(current_gpl_info$series_id)
length(current_gpl_info$sample_id)

# Accessing expression data, in the supplemental files

sfiles = getGEOSuppFiles("GSE158890")
fnames = rownames(sfiles)

data = read.delim(fnames[1], header=TRUE, check.names = FALSE)

summarized_gene_counts <- sort(table(data$gene_name), decreasing = TRUE)

summarized_gene_counts[1:7]


cpms = cpm(data[,9:length(colnames(data))])
rownames(cpms) <- data[,1]
keep = rowSums(cpms >1) >=10
data_filtered = data[keep,]

summarized_gene_counts_filtered <- sort(table(data_filtered$gene_name), decreasing = TRUE)
summarized_gene_counts_filtered[ which(summarized_gene_counts_filtered>1)[1:10]]





#normie

data2plot <- log2(cpm(data_filtered[,9:length(colnames(data_filtered))]))
boxplot(data2plot, xlab = "Samples", ylab = "log2 CPM",
        las = 2, cex = 0.5, cex.lab = 0.5,
        cex.axis = 0.5, main = "HIF/NPM1 RNASeq Samples")
#draw the median on each box plot
abline(h = median(apply(data2plot, 2, median)),
       col = "green", lwd = 0.6, lty = "dashed")





samples <- data.frame(lapply(colnames(data)[,9:length(colnames(data))],
                             FUN=function(x){unlist(strsplit(x, split = "\\."))[c(2,3)]}))
colnames(samples) <- colnames(data)[,9:length(colnames(data))]
rownames(samples) <- c("cohort","cell_type")
samples <- data.frame(t(samples))


filtered_data_matrix <- as.matrix(data_filtered[,9:length(colnames(data_filtered))])
rownames(filtered_data_matrix) <- data_filtered$gene_id
d = DGEList(counts=filtered_data_matrix, group=samples$cell_type)


d = calcNormFactors(d)


























#Converting HUGOs

convertedHugos <- "hugos.rds"

if(file.exists(convertedHugos)) {

  converted <- readRDS(convertedHugos)

} else {

  mart <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

  genes <- c(data_filtered$gene_id)

  converted <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes, mart = mart)

  hugos <- converted$hgnc_symbol

  saveRDS(converted, convertedHugos)

}

#Merginging successfully converted HUGOs with original Ensembl list. Some will not have matched.

x <- merge(converted, data_filtered, by.x = 1, by.y = 0, all.y=TRUE)

missing <- x$ensembl_gene_id[ which(is.na(x$hgnc_symbol))]

length(missing)

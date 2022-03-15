if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

if (!requireNamespace("GEOmetadb", quietly = TRUE))
  BiocManager::install("GEOmetadb", force = TRUE)
library(GEOmetadb)

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")
library('GEOquery')

if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")ge
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

y = getGEOSuppFiles("GSE158890")

fnames = rownames(y)

hif_exp = read.delim(fnames[1], header=TRUE, check.names = FALSE)

hif_exp[1:15, 4:17]

samples <- data.frame(lapply(colnames(hif_exp)[9:18],
                             FUN=function(x){unlist(strsplit(x,
                                                             split = "\\_"))[c(1,3)]}))
colnames(samples) <- colnames(hif_exp)[9:18]
rownames(samples) <- c("cell_type","patients")
samples <- data.frame(t(samples))


summarized_gene_counts <- sort(table(hif_exp$gene_name), decreasing = TRUE)

summarized_gene_counts[1:7]


cpms = cpm(hif_exp[,9:length(colnames(hif_exp))])
rownames(cpms) <- hif_exp[,1]
keep = rowSums(cpms >1) >=10
hif_exp_filtered = hif_exp[keep,]

summarized_gene_counts_filtered <- sort(table(hif_exp_filtered$gene_name), decreasing = TRUE)

summarized_gene_counts_filtered[ which(summarized_gene_counts_filtered>1)[1:10]]


filtered_data_matrix <- as.matrix(hif_exp_filtered[,9:18])
rownames(filtered_data_matrix) <- hif_exp_filtered$gene_id

d = DGEList(counts=filtered_data_matrix, group=samples$cell_type)
d = calcNormFactors(d)
normalized_counts <- cpm(d)

plotMDS(d, labels=rownames(samples),
        col = c("darkgreen","blue")[factor(samples$cell_type)])


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

x <- x$ensembl_gene_id[ which(is.na(x$hgnc_symbol))]

length(missing)


# IDK





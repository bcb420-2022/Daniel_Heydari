# GSEA Homework BCB420Y

#install required R and bioconductor packages
tryCatch(expr = { library("RCurl")},
         error = function(e) {  install.packages("RCurl")},
         finally = library("RCurl"))

tryCatch(expr = { library("BiocManager")},
         error = function(e) {
           install.packages("BiocManager")},
         finally = library("BiocManager"))
tryCatch(expr = { library("ggplot2")},
         error = function(e) { install.packages("ggplot2")},
         finally = library("ggplot2"))
#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library("RCy3")},
         error = function(e) { BiocManager::install("RCy3")},
         finally = library("RCy3"))

#path to GSEA jar
# In order to run GSEA automatically you need to speciry the path to the gsea jar file.
#With the latest release of gsea (4.0.2) they no longer release a bundled jar
# and instead release a scriptted way to launch the gsea client.
# specify the java version as 11 if you are using the later version gsea
# the gsea_jar also needs to be the full path to the GSEA 4.0.2 directory that you
# downloaded from GSEA. for example (/Users/johnsmith/GSEA_4.0.2/gsea-cli.sh)
gsea_jar <- params$gsea_jar
java_version <- params$java_version
#Gsea takes a long time to run.  If you have already run GSEA manually or previously there is no need to re-run GSEA.  Make sure the
# gsea results are in the current directory and the notebook will be able to find them and use them.
run_gsea = params$run_gsea
#navigate to the directory where you put the downloaded protocol files.
working_dir <- params$working_dir
# leave blank if you want the notebook to discover the gsea directory for itself
#gsea_directory = paste(working_dir,"Mesen_vs_Immuno.GseaPreranked.1497635459262",sep="/")
gsea_directory = params$gsea_directory
analysis_name <- params$analysis_name
rnk_file <- params$rnk_file
expression_file <- params$expression_file
classes_file <- params$classes_file

gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"
# list all the files on the server
filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)
# get the gmt that has all the pathways and does not include terms inferred
# from electronic annotations(IEA) start with gmt file that has pathways only
rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)", contents,
              perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))
dest_gmt_file <- file.path(getwd(), gmt_file)
download.file(paste(gmt_url, gmt_file, sep = ""), destfile = dest_gmt_file)

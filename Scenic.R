library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)

###input gene exrpession matrix###
library(SingleCellExperiment)
sheepRumen <- read.table("SCENIC/exprMat.txt")
exprMat<-as.matrix(sheepRumen)
dim(exprMat)
mydbDIR="/SCENIC/Resource/hg38_scenic/"
dir(mydbDIR)
mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org = "hgnc", 
                                  nCores = 4,
                                  dbDir = mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNC")
saveRDS(scenicOptions, "SCENIC/int/scenicOptions.rds")

###gene filtering##
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 0.015 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)

###Calculating correlation matrix###
runCorrelation(exprMat_filtered, scenicOptions)

###Calculating TF-targets correlation###
runGenie3(exprMat_filtered, scenicOptions, nParts = 20)
runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions@settings$nCores <- 2

###Inferring transcriptional regulatory networks##
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat)
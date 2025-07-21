#!/usr/bin/Rscript

#######====================== MPRAnalyze_mgx.R ======================#######
#MPRAnalyze_mgx.R is a wrapper script for applying MPRAnalyze to paired 
#MGX and MTX data using feature level MGX counts as the DNA counts. 

#MPRAnalyze Model specifications: 

#mgx: dna.counts used for the DNA model are feature-level metagenomic counts

#Depth factors are estimated within each sample assuming that each MGX/MTX 
#sample represents an individual library. This is done using 
#lib.factor = "SampleID" which must be a column in metadata containing
#unique values for each sample.  

#Model design formula only includes experimental condition ("Phenotype") 
#for both the dna and rna models; barcodes are assumed to be one per sample
#(i.e. no technical replicates of any library)

#Arguments: 
# dataset_name : file handle for input files in input_dir. Input files should be of format
#                {dataset_name}.{bug/mgx/mtx}_abunds.tsv etc.  
# input_dir: file path to directory containing input files (naming convention described for dataset_name).
#             The paths {input_dir}/{dataset_name}.{bug/mgx/mtx}_abunds.tsv must exist 
# output_dir: file path to directory which will be created for DGE results files 

#The below R code was modified from the MPRAnalyze vignette:
#https://bioconductor.org/packages/release/bioc/vignettes/MPRAnalyze/inst/doc/vignette.html

#The synthetic data can be found under Simulated Datasets: http://huttenhower.sph.harvard.edu/mtx2021

#Last Modified: 10/10/23
#CHANGE LOG 
#10/10/24: Script created based on MPRAnalyze_test and MPRAnalyze_synth_test.
#######================================================================#######

#Set dataset_name and input_dir, for loading input data and 
#output_dir, where DGE results will be saved 
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<=2) {
  stop("Please provide three arguments: dataset_name, input_dir, and output_dir.")
} else {
  dataset_name <- args[1]
  input_dir <- args[2]
  output_dir <- args[3] 
}

library(MPRAnalyze)
#File paths for counts and metadata inputs
dna.fpath <-file.path(input_dir,sprintf("%s.mgx_abunds.tsv",dataset_name)) #mgx_abunds (feature-level) for M6
rna.fpath <- file.path(input_dir,sprintf("%s.mtx_abunds.tsv",dataset_name))
metadata.fpath <-file.path(input_dir,sprintf("%s.MPRA_metadata.tsv",dataset_name))
#Read in counts data and metadata 
dna.counts <- read.table(dna.fpath)
rna.counts <- read.table(rna.fpath)
metadata <- read.table(metadata.fpath)

#do not provide controls (vector of bool labels for control enhancers)
obj <- MpraObject(dnaCounts = as.matrix(dna.counts) , rnaCounts = as.matrix(rna.counts), 
                  dnaAnnot = metadata, rnaAnnot = metadata)
#unique values of "SampleID" specify individual dna/rna libraries (samples)
#across which depth normalization is performed
#use default upper-quartile normalization factors
obj <- estimateDepthFactors(obj, lib.factor = c("SampleID"),
                            which.lib = "both")
#dnaDesign and rnaDesign only contain one variable, Phenotype
#reducedDesign represents the null hypothesis (no difference between groups)
obj <- analyzeComparative(obj = obj, 
                          dnaDesign = ~ Phenotype, 
                          rnaDesign = ~ Phenotype, 
                          reducedDesign = ~ 1)
#Use likelihood ratio test for statistical testing for DE 
res <- testLrt(obj)
#Order results by fdr and write to file 
res_ordered <- res[order(res$fdr),]
res_fpath = file.path(output_dir,"all_results.tsv")
write.table(res_ordered,res_fpath,quote=FALSE,sep='\t')
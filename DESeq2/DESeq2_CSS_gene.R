#!/usr/bin/Rscript

#######====================== DESeq2_CSS_gene.R ======================#######
#DESeq2_CSS_gene uses DESeq2 with community level scaling and gene-level DNA covariates
#(provided as a normalization matrix instead of as part of the design matrix). 

#arguments
# dataset_name : file handle for input files in input_dir. Input files should
#   be of format {dataset_name}.{bug/mgx/mtx}_abunds.tsv etc.
# input_dir: file path to directory containing input files (naming convention described for dataset_name)
# output_dir: file path to directory which will be created for DGE results files 

#Author: Evan Lee
#Last updated: 04/30/25
# CHANGE LOG
# 04/25: implemented minimum dispersion error handling. 
#   Also fixed bug in normFactor calculation. Need to redo synth results 
#   with updated script. 
#######================================================================#######

#Load packages 
suppressPackageStartupMessages(require('DESeq2', character.only = TRUE))
suppressPackageStartupMessages(require('dplyr', character.only = TRUE))
library(stats)
# library(apeglm) #for lfcShrink using apeglm prior - opting for ashr because allows for 
# use of contrasts rather than needing explicit factor level setting in col_data.
library(ashr) #for lfcShrink using ashr prior
#Set dataset_name and input_dir for loading input data and output_dir where DGE results will be saved 
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<=2) {
  stop("Please provide three arguments: dataset_name, input_dir, and output_dir.")
} else {
  dataset_name <- args[1]
  input_dir <- args[2]
  output_dir <- args[3] 
}
dir.create(output_dir,showWarnings = FALSE,recursive = TRUE) #make output_dir if doesn't exist yet 

#File paths 
bug_fpath <- file.path(input_dir,sprintf("%s.bug_abunds.tsv",dataset_name))
mgx_fpath <- file.path(input_dir,sprintf("%s.mgx_abunds.tsv",dataset_name))
mtx_fpath <- file.path(input_dir,sprintf("%s.mtx_abunds.tsv",dataset_name))
md_fpath <- file.path(input_dir,sprintf("%s.metadata.tsv",dataset_name))

#Read in input data 
bug.df <- read.table(bug_fpath,header = TRUE,row.names=1)
mgx.df <- read.table(mgx_fpath,header=TRUE,row.names=1)
mtx.df <- read.table(mtx_fpath,header=TRUE,row.names=1)
metadata.df <- read.table(md_fpath,header=TRUE,row.names=1)

#condition covariate data, formatted for DESeq2. 
col_data <- as.data.frame(t(metadata.df['Phenotype',]))
#Set factor levels such that phenotype 0 (glucose) is reference
# This is important because lfcShrink using apeglm priors does not allow specification of contrast
col_data$Phenotype <- factor(col_data$Phenotype,levels=c(0,1))

#Make DESeq dataset for entire mtx.df with just col_data and Phenotype design
all_dds <- DESeqDataSetFromMatrix(mtx.df,colData=col_data,design=~ Phenotype)
##Make norm_matrix with gene level normalization factors (metagenomic DNA abundances of features)
#Gene-level normalization factors: https://support.bioconductor.org/p/90817/
PSEUDOCOUNT <- 1
norm_matrix <- (mgx.df + PSEUDOCOUNT) / exp(rowMeans(log(mgx.df + PSEUDOCOUNT)))
#estimate size factors for this dds for community-level scaling
#Note still need to use poscounts size factor estimation since all genes have at least one zero 
all_dds <- estimateSizeFactors(all_dds, normMatrix=norm_matrix, type = 'poscounts')

#In CSS-none, do not need to iterate by each bug; fit dispersions and log fold change trends 
#over entire dataset as might be done in a naive MTX approach  
#Run DGE analysis, implementing minimum dispersion error handling. 

all_dds <- tryCatch({DESeq(all_dds)},
  error = function(e){
    #Brief: default estimateDispersions within DESeq will fail to fit a trend
    # when all genes from a taxon/bug are very low expression (i.e. for low
    # abundance organisms/MAGs) 
    #Message: "all gene-wise dispersion estimates are within 2 orders of 
    # magnitude from the minimum value, and so the standard curve fitting 
    # techniques will not work."
    #Developer recommendations: use estimateDispersionsGeneEst (i.e. 
    # no EB dispersion shrinkage), set dispersions for dds, and then 
    # perform statistical testing with preferred test (nbinomWaldTest)
    if (grepl("all gene-wise dispersion estimates are within 2 orders of magnitude", e$message)){
      #If the error is related to dispersion fitting, use gene-wise estimates
      # as final estimates as recommended by error message. 
      all_dds <- estimateDispersionsGeneEst(all_dds)
      dispersions(all_dds) <- mcols(all_dds)$dispGeneEst
      all_dds <- nbinomWaldTest(all_dds)
      return(all_dds)
    }
  }
)
#Formatting results for ashr logFC shrinkage using contrasts:
#Apply logFC shrinkage to results; using ashr over apeglm 
#because allows for use of contrasts and methods are otherwise
#equivalent/recommended by developers for shrinkage.
# Not providing res; with type='ashr', res overrides contrast 
all_res <- lfcShrink(all_dds,
                      contrast=c('Phenotype','1','0'),
                      type="ashr")

#Sort by padj and write output 
all_res$padj_all <- stats::p.adjust(all_res$pvalue,method="BH")
all_res_ordered <- all_res[order(all_res$padj_all),]
all_res_fpath <- file.path(output_dir,"all_results.tsv")
write.table(all_res_ordered,all_res_fpath,quote = FALSE,sep='\t')

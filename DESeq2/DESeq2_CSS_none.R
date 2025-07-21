#!/usr/bin/Rscript

#######====================== DESeq2_CSS_none.R ======================#######
#DESeq2_CSS_none uses DESeq2 with community level scaling and no DNA covariates. 

#arguments
# dataset_name : file handle for input files in input_dir. Input files should be of format
#                {dataset_name}.{bug/mgx/mtx}_abunds.tsv etc.  
# input_dir: file path to directory containing input files (naming convention described for dataset_name)
# output_dir: file path to directory which will be created for DGE results files 

#Author: Evan Lee
#Last updated: 8/20/24
# CHANGE LOG
# 05/25: implemented minimum dispersion error handling, updated for consistency
#         with TSS  version. 
# 8/20/24: updated for BG01 - removed variables for mgx/bug abunds
#           Desired call is Rscript --vanilla DESeq2_CSS_none.R pco100 . deseq2_output 
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
mtx_fpath <- file.path(input_dir,sprintf("%s.mtx_abunds.tsv",dataset_name))
md_fpath <- file.path(input_dir,sprintf("%s.metadata.tsv",dataset_name))

#Read in input data 
mtx.df <- read.table(mtx_fpath,header=TRUE,row.names=1)
metadata.df <- read.table(md_fpath,header=TRUE,row.names=1)

#condition covariate data, formatted for DESeq2. 
col_data <- as.data.frame(t(metadata.df['Phenotype',]))
#Set factor levels such that phenotype 0 (glucose) is reference
# This is important because lfcShrink using apeglm priors does not allow specification of contrast 
col_data$Phenotype <- factor(col_data$Phenotype,levels=c(0,1))

#Make DESeq dataset for entire mtx.df with just col_data; 
all_dds <- DESeqDataSetFromMatrix(mtx.df,colData=col_data,design=~ Phenotype)
#estimate size factors for this dds for community-level scaling
#Note still need to use poscounts size factor estimation since all genes have at least one zero 
all_dds <- estimateSizeFactors(all_dds, type = 'poscounts')
all_size_factors <- sizeFactors(all_dds) #equivalent to dds@colData@listData$sizeFactor used for TSS Klingenberg and Meinicke

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
      # If the error is related to dispersion fitting, use gene-wise estimates as final estimates as recommended by error message
      all_dds <- estimateDispersionsGeneEst(all_dds)
      dispersions(all_dds) <- mcols(all_dds)$dispGeneEst
      all_dds <- nbinomWaldTest(all_dds)
      return(all_dds)
    }
  }
)
# all_res <- results(all_dds,contrast=c('Phenotype','1','0'))
#Apply logFC shrinkage; using ashr because allows specification of contrast 
all_res <- lfcShrink(all_dds,
                      contrast=c('Phenotype','1','0'),
                      type="ashr")

#Sort by padj and write output 
all_res$padj_all <- stats::p.adjust(all_res$pvalue,method="BH")
all_res_ordered <- all_res[order(all_res$padj_all),]
all_res_fpath <- file.path(output_dir,"all_results.tsv")
write.table(all_res_ordered,all_res_fpath,quote = FALSE,sep='\t')

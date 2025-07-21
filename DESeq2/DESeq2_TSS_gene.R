#!/usr/bin/Rscript

#######====================== DESeq2_TSS_gene.R ======================#######
#DESeq2_TSS_gene uses DESeq with taxon-specific scaling and gene-level DNA covariates
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
# 9/24: updated for BG01
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

#PSEUDOCOUNT is used in generating the gene-level normalization factors matrix (norm_matrix)
PSEUDOCOUNT <- 1

#condition covariate data, formatted for DESeq2. 
col_data <- as.data.frame(t(metadata.df['Phenotype',]))
#Set factor levels such that phenotype 0 (glucose) is reference
# This is important because lfcShrink using apeglm priors does not allow specification of contrast
col_data$Phenotype <- factor(col_data$Phenotype,levels=c(0,1))
#Table for all statistical test results 
all_res <- data.frame()

#Iterate through each bug to apply taxon-specific scaling
for (bug in row.names(bug.df)) {
  #Select features for this bug with filter, grepl
  bug_mgx <- mgx.df %>% filter(grepl(bug,row.names(mgx.df)))
  bug_mtx <- mtx.df %>% filter(grepl(bug,row.names(mtx.df))) 
  #A separate bug_col_data variable is used however to facilitate zero-filtering of samples 
  #with zero counts for all of the bug genes
  bug_col_data <- col_data #In gene-covariate TSS, the design matrix doesn't change per bug (only contains condition/phenotype information)
  #https://support.bioconductor.org/p/101179/ = check that column sums all non-zero
  #Causes na values in size factor estimation hence filtering - remove samples 
  #with all zero counts for genes from this bug
  #Filter out columns (samples) from mgx, mtx, and col data which
  #contain all zeros in mtx data
  if (0 %in% colSums(bug_mtx)) {
    zero_counts_samples <- colnames(bug_mtx)[colSums(bug_mtx)==0]
    bug_mgx <- bug_mgx %>% select(!all_of(zero_counts_samples))
    bug_mtx <- bug_mtx %>% select(!all_of(zero_counts_samples))
    bug_col_data <- col_data %>% filter(!row.names(col_data) %in% zero_counts_samples)
    
    #Check whether zero-filtered bug_col_data Phenotype still contains samples from both groups 
    if (length(unique(bug_col_data$Phenotype))==1) {
      print(sprintf("Filtering out samples with zeros in all rows only leaves one experimental condition for bug %s.", bug))
      print(sprintf("Skipping DGE testing for bug %s.",bug))
      next
    } 
  }
  #Generate DESeqDataSet 
  bug_dds <- DESeqDataSetFromMatrix(bug_mtx,
                                    colData = bug_col_data,
                                    design= ~ Phenotype)
  #Gene-level normalization factors: https://support.bioconductor.org/p/90817/
  # norm_matrix <- bug_mgx #Fails because 0s in bug_mgx and resulting counts/normMatrix has na
  #Add pseudocount and normalize by geometric mean as recommended by Love for normMatrix 
  norm_matrix <- (bug_mgx + PSEUDOCOUNT) / exp(rowMeans(log(bug_mgx + PSEUDOCOUNT)))
  #poscounts sizeFactor estimation tolerates if every feature has a 0 
  #Gene-level normalization factors: https://support.bioconductor.org/p/90817/
  bug_dds <- estimateSizeFactors(bug_dds, normMatrix=norm_matrix,type = 'poscounts')
  
  #DGE Testing - Implementing Marie's handling of dispersion fitting errors.
  # If no error, tryCatch will run DESeq(bug_dds) as normal. 
  bug_dds <- tryCatch({DESeq(bug_dds)},
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
        bug_dds <- estimateDispersionsGeneEst(bug_dds)
        dispersions(bug_dds) <- mcols(bug_dds)$dispGeneEst
        bug_dds <- nbinomWaldTest(bug_dds)
        return(bug_dds)
      }
    }
  )
  #Formatting results for ashr logFC shrinkage using contrasts:
  #Apply logFC shrinkage to results; using ashr over apeglm 
  #because allows for use of contrasts and methods are otherwise
  #equivalent/recommended by developers for shrinkage.
  # Not providing res; with type='ashr', res overrides contrast 
  res <- lfcShrink(bug_dds,
                      contrast=c('Phenotype','1','0'),
                      type="ashr")
  #Concatenate bug results to all results
  all_res <- rbind(all_res,data.frame(res))
}
#Sort by padj and write output 
all_res$padj_all <- stats::p.adjust(all_res$pvalue,method="BH")
all_res_ordered <- all_res[order(all_res$padj_all),]
all_res_fpath <- file.path(output_dir,"all_results.tsv")
write.table(all_res_ordered,all_res_fpath,quote = FALSE,sep='\t')


#!/usr/bin/Rscript

#######====================== M5_lenient_MTX.R ======================#######
#M5_lenient_MTX is a wrapper script for MTX_model using taxon-level DNA covariates. Lenient (feature-level) zero filtering is used. 

#Model specifcations: log transform, total sum scaling (TSS) normalization, fixed_effects of Phenotype (in metadata file).
#Model specifcations (continued): using per-taxon DNA abundances in input_dnadata. 'lenient' rna_dna_filter 

#Arguments: 
# dataset_name : file handle for input files in input_dir. Input files should be of format
#                {dataset_name}.{bug/mgx/mtx}_abunds.tsv etc.  
# input_dir: file path to directory containing input files (naming convention described for dataset_name).
#             The paths {input_dir}/{dataset_name}.{bug/mgx/mtx}_abunds.tsv must exist 
# output_dir: file path to directory which will be created for DGE results files 

#The below R code was based on https://github.com/biobakery/MTX_model/blob/master/README.md
#The synthetic data can be found under Simulated Datasets: http://huttenhower.sph.harvard.edu/mtx2021

#note that the original file (synth_mgx_mtx.tar.gz) should be unpacked with tar xvf synth_mgx_mtx.tar.gz 
#Last Modified: 9/26/23
#CHANGE LOG 
#09/26/24: added additional optioons, most importantly to override prevalence filtering and try to align to Zhang published results 
#in terms of the number of features tested 
#09/11/23: split into strict vs lenient scripts with different 'rna_dna_filt' selections
#######================================================================#######

#Set dataset_name and input_dir for loading input data and output_dir where DGE results will be saved 
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<=2) {
  stop("Please provide three arguments: dataset_name, input_dir, and output_dir.")
} else {
  dataset_name <- args[1]
  input_dir <- args[2]
  output_dir <- args[3] 
}

library(MTXmodel)
input_data <- file.path(input_dir,sprintf("%s.mtx_abunds.tsv",dataset_name))
input_metadata <-file.path(input_dir,sprintf("%s.metadata.tsv",dataset_name))
input_dnadata <-file.path(input_dir,sprintf("%s.mgx_bug_abunds.tsv",dataset_name)) #bug_abunds for M5
fit_data <- MTXmodel(
    input_data, input_metadata, output_dir, transform = "LOG",
    min_prevalence = 0, #override default min_prevalence of 0.1 which seems incongruent with what is used for their published test results
    max_significance = 0.05, #override default alpha = 0.25; only affects significance filtering in significant_results.tsv
    fixed_effects = c('Phenotype'),
    reference = 0,
    normalization = 'TSS',
    standardize = FALSE,
    input_dnadata = input_dnadata,
    rna_dna_flt = "lenient",
    plot_heatmap = FALSE,
    plot_scatter = FALSE #skip scatter/boxplots from MaAsLin
    )
#!/usr/bin/env Rscript

setwd('/Users/leizhang/Desktop/shiLab/BreastCancer')

## Split data
gpr_dir = './gpr'
map_file = './block_ind_mappint.csv'	
out_dir = './originalData'
ID_to_change = NULL # NULL if no protein is subjected to the id change 
source('bin/splitData.R')

## Preprocess
input_dir = out_dir
delNames = c("buffer", "cy3/cy5", "Human IgG", "Human IgM")
delSampleNames = NULL # NULL if no bad samples are deleted
delFlags_minus50 = TRUE
delFlags_minus100 = FALSE 
outfile = './raw_allSamples.csv'
normalized_outfile =  './normalization/normalized_allSamples.csv'
source('bin/preprocess.R')

## Test
source('bin/difftest.R')
group_inputfile = './GroupInfo.csv'
cmpInfo = c('ZvsA', 'AvsB', 'BvsC', 'ZvsA_B_C', 'AvsB_C', 'ZvsB_C')
logTransform = T
outputPath = './comparison'
rel_ttest = diffTest(norm_inputfile=normalized_outfile, group_inputfile, cmpInfo, logTransform=T, outputPath=outputPath)

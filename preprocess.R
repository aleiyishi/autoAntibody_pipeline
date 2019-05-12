#!/usr/bin/env Rscript

# ChangeLog:
## April 18, 2019: add the option to dealing with the flag column

library(ggplot2)
library(reshape2)
library(pheatmap)
library(limma)

## extract the specified columns and calculate the F635 value with background subtracted

extract_cols <- function(listObj, cols=c("ID", "F635 Median", "B635 Median", 'Flags')){
	
	out_list = list()

	for (m in 1:length(listObj)){
		gpr = listObj[[m]][cols]
		gpr$F635_subBg = NA
		F635_subBg = gpr[,cols[2]] - gpr[,cols[3]]
		F635_subBg2 = ifelse(F635_subBg<0, 0, F635_subBg)		
		
		if (delFlags_minus50 & !delFlags_minus100) {
			
			message(sprintf('[Note]: Intensity with flags -50 will be set to NA in %s', names(listObj)[m]))
			F635_subBg3 = ifelse(as.character(gpr[,cols[4]]) == '-50', NA, F635_subBg2)
			gpr$F635_subBg = F635_subBg3
		
		} else if (delFlags_minus50 & delFlags_minus100) {

			message(sprintf('[Note]: Intensity with flags -50 or -100 will be set to NA in %s', names(listObj)[m]))
			F635_subBg3 = ifelse(as.character(gpr[,cols[4]]) == '-50' | as.character(gpr[,cols[4]]) == '-100' , NA, F635_subBg2)
			gpr$F635_subBg = F635_subBg3

		}

		gpr_split = split(gpr[c(cols[1], 'F635_subBg')], f=gpr[,cols[1]])
		gpr_numeric = sapply(gpr_split, function(x)mean(x[,2], na.rm=T))
		out_list[[names(listObj)[m]]] = gpr_numeric
			
	}
	
	return(do.call(cbind, out_list))

}

gpr_files = list.files(input_dir, full.names=T)
gpr_data_list = lapply(gpr_files, function(x)read.table(x, header=T, sep='\t', check.names=F, as.is=T))
names(gpr_data_list) = sub('.txt', '', basename(gpr_files))
gprData = extract_cols(gpr_data_list)
gprData = gprData[!rownames(gprData) %in% delNames,!colnames(gprData) %in% delSampleNames]
write.csv(gprData, outfile, quote=F)

## perform the quantile normalization
normalized_dir = dirname(normalized_outfile)
if (!dir.exists(normalized_dir)) dir.create(normalized_dir)

## raw data
gprData_long = melt(gprData)
pdf(file.path(normalized_dir, 'boxplot_signal_distribution_before_normalization.pdf'), width=42, height=12)
gBefore = ggplot(gprData_long, aes(x=as.factor(Var2),y=value))+geom_boxplot(outlier.size=0.5, outlier.shape=NA)+labs(x='',y='Raw expression level') + theme(axis.text.x=element_text(size=6, angle=90))
print(gBefore)
dev.off()

## normalized data
gprData_normalized = normalizeBetweenArrays(gprData, method='quantile')
write.csv(gprData_normalized, file=normalized_outfile, quote=F)
write.csv(log2(gprData_normalized), file=sub('.csv', '_logTransformed.csv', normalized_outfile), quote=F)

gprData_normalized_long = melt(gprData_normalized)
pdf(file.path(normalized_dir, 'boxplot_signal_distribution_after_quantile_normalization.pdf'), width=42, height=12)
gAfter = ggplot(gprData_normalized_long, aes(x=as.factor(Var2),y=value))+geom_boxplot(outlier.size=0.5, outlier.shape=NA)+labs(x='',y='Normalized expression level') + theme(axis.text.x=element_text(size=6, angle=90))
print(gAfter)
dev.off()

## draw the heatmap
### raw data
pdf('HC_clustering_of_all_samples_based_on_raw_data.pdf', width=20)
p <- pheatmap(gprData, annotation_col = NA, show_rownames=F, scale='row', fontsize_col=3, cutree_cols=2)
dev.off()




















#!/usr/bin/env Rscript
## split the log2-normalized data by the group such as A, B, C and so on

library(pheatmap)

norm_inputfile =  './normalization/normalized_allSamples_logTransformed.csv'
group_inputfile = './GroupInfo.csv'
group_tobeClustered = c('Z', 'A', 'B', 'C')

#######################################################################################
normData = read.csv(norm_inputfile, as.is=T, check.names=F, row.names=1)
groupInfo = read.csv(group_inputfile, as.is=T, check.names=F)

groupList = split(groupInfo$ID, f=groupInfo$Group)

for (k in 1:length(groupList)){

	groupName = names(groupList)[k]
	groupData = normData[,colnames(normData) %in% groupList[[k]]]
	write.csv(groupData, sub('.csv', sprintf('_%s.csv', groupName), norm_inputfile), quote=F)

}



### log2 transformed
#unionPro = read.table('bin/ZvsBC.txt', as.is=T)
#group_tobeClustered = c('Z', 'C')
sample_tobeClustered = groupInfo$ID[groupInfo$Group %in% group_tobeClustered]
normData_tobeClustered = normData[,colnames(normData) %in% sample_tobeClustered]
annotation_col = groupInfo[groupInfo$ID %in% colnames(normData_tobeClustered),]
annotation_col_2 = annotation_col['Group']
rownames(annotation_col_2) = annotation_col$ID
#pdf('HC_clustering_of_all_samples_based_on_log2-normalized_data.pdf', width=20)
p <- pheatmap(normData_tobeClustered[unionPro[,1],], annotation_col = annotation_col_2, show_rownames=F, scale='none', fontsize_col=3,  clustering_distance_cols = "correlation", clustering_method = "average")
#dev.off()
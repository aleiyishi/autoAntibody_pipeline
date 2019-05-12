#!/usr/bin/env Rscript
## split the log2-normalized or raw data by the group such as A, B, C and so on

#library(pheatmap)
args = commandArgs(TRUE)
norm_inputfile = args[1]
group_inputfile = args[2]

## test
#norm_inputfile =  './normalization/normalized_allSamples_logTransformed.csv'
#group_inputfile = './GroupInfo.csv'

#######################################################################################
normData = read.csv(norm_inputfile, as.is=T, check.names=F, row.names=1)
groupInfo = read.csv(group_inputfile, as.is=T, check.names=F)

groupList = split(groupInfo$ID, f=groupInfo$Group)

for (k in 1:length(groupList)){

	groupName = names(groupList)[k]
	groupData = normData[,colnames(normData) %in% groupList[[k]]]
	write.csv(groupData, sub('.csv', sprintf('_%s.csv', groupName), norm_inputfile), quote=F)

}

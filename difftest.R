#Released in April 3th, 2019

diffTest <- function(norm_inputfile, group_inputfile, cmpInfo, logTransform=T, outputPath='.'){

	ttestObj = list()

	normData = read.csv(norm_inputfile, as.is=T, check.names=F, row.names=1)
	groupInfo = read.csv(group_inputfile, as.is=T, check.names=F)
	cmnID = intersect(colnames(normData), groupInfo$ID)
	groupInfo_cmn = groupInfo[groupInfo$ID %in% cmnID, ]

	for (i in 1:length(cmpInfo)){
		
		path = file.path(outputPath, cmpInfo[i])
		dir.create(path, recursive=T)
		message(sprintf('Perform the T test bewteen %s', cmpInfo[i]))
		cmpGroup = unlist(strsplit(cmpInfo[i], 'vs'))
		caseLab = unlist(strsplit(cmpGroup[2], '_'))
		ctrLab = unlist(strsplit(cmpGroup[1], '_'))
		caseID = groupInfo_cmn$ID[groupInfo_cmn$Group %in% caseLab]
		ctrID = groupInfo_cmn$ID[groupInfo_cmn$Group %in% ctrLab]
		caseData = normData[,caseID]
		ctrData = normData[,ctrID]

		tin_list = list()
		for (k in 1:nrow(normData)){

			if (!logTransform) {
				rel_ttest = t.test(ctrData[k,], caseData[k,])
				simpleRel = cbind(Name=rownames(normData)[k], Pvalue=rel_ttest$p.value, Mean_control=unname(rel_ttest$estimate[1]), Mean_case=unname(rel_ttest$estimate[2]), FC=as.numeric(unname(rel_ttest$estimate[2]))/as.numeric(unname(rel_ttest$estimate[1])))
				tin_list[[k]] = simpleRel
			} else {
				rel_ttest = t.test(log2(ctrData[k,]), log2(caseData[k,]))
				simpleRel = cbind(Name=rownames(normData)[k], Pvalue=rel_ttest$p.value, Mean_control=unname(rel_ttest$estimate[1]), Mean_case=unname(rel_ttest$estimate[2]), FC=2^(as.numeric(unname(rel_ttest$estimate[2]))-as.numeric(unname(rel_ttest$estimate[1]))))
				tin_list[[k]] = simpleRel
			}
			
		}
		tin = do.call(rbind, tin_list)
		write.csv(tin, file=file.path(path, 'Result_ttest.csv'), quote=F, row.names=F)
		ttestObj[[i]] = tin


		## cluster sampls based on differential markers
		if (logTransform) {

			p1 = 0.05
			diffMarkers = tin[,'Name'][as.numeric(tin[,'Pvalue']) <= p1]

			if (length(diffMarkers) > 2) {
				normData_tobeclustered = cbind(log2(ctrData), log2(caseData))[diffMarkers,]
				annotation_col = data.frame(Group=factor(c(rep(cmpGroup[1], ncol(ctrData)), rep(gsub('_', '+', cmpGroup[2], fixed=T), ncol(caseData)))), stringsAsFactors=F)
				rownames(annotation_col) = colnames(normData_tobeclustered)
				pdf(file.path(path, sprintf('HC_clustering_of_%s_samples_based_on_log2-normalized_data_%s.pdf', cmpInfo[i], as.character(p1))), width=nrow(annotation_col)/10)
				hc_p1 <- pheatmap::pheatmap(normData_tobeclustered, annotation_col = annotation_col, show_rownames=T, scale='none', fontsize_col=4.5, clustering_distance_cols = "euclidean", clustering_method = "average")
				print(hc_p1)
				dev.off()
			}
			
			p2 = 0.01
			diffMarkers = tin[,'Name'][as.numeric(tin[,'Pvalue']) <= p2]

			if (length(diffMarkers) > 2) {
				normData_tobeclustered = cbind(log2(ctrData), log2(caseData))[diffMarkers,]
				annotation_col = data.frame(Group=factor(c(rep(cmpGroup[1], ncol(ctrData)), rep(gsub('_', '+', cmpGroup[2], fixed=T), ncol(caseData)))), stringsAsFactors=F)
				rownames(annotation_col) = colnames(normData_tobeclustered)
				pdf(file.path(path, sprintf('HC_clustering_of_%s_samples_based_on_log2-normalized_data_%s.pdf', cmpInfo[i], as.character(p2))), width=nrow(annotation_col)/10)
				hc_p2 <- pheatmap::pheatmap(normData_tobeclustered, annotation_col = annotation_col, show_rownames=T, scale='none', fontsize_col=4.5, clustering_distance_cols = "euclidean", clustering_method = "average")
				print(hc_p2)
				dev.off()
			}

			p3 = 0.001
			diffMarkers = tin[,'Name'][as.numeric(tin[,'Pvalue']) <= p3]

			if (length(diffMarkers) > 2) {
				normData_tobeclustered = cbind(log2(ctrData), log2(caseData))[diffMarkers,]
				annotation_col = data.frame(Group=factor(c(rep(cmpGroup[1], ncol(ctrData)), rep(gsub('_', '+', cmpGroup[2], fixed=T), ncol(caseData)))), stringsAsFactors=F)
				rownames(annotation_col) = colnames(normData_tobeclustered)
				pdf(file.path(path, sprintf('HC_clustering_of_%s_samples_based_on_log2-normalized_data_%s.pdf', cmpInfo[i], as.character(p3))), width=nrow(annotation_col)/10)
				hc_p3 <- pheatmap::pheatmap(normData_tobeclustered, annotation_col = annotation_col, show_rownames=T, scale='none', fontsize_col=4.5, clustering_distance_cols = "euclidean", clustering_method = "average")
				print(hc_p3)
				dev.off()
			}


		}

	}
	
	return(ttestObj)

}




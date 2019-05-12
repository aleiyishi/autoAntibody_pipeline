#!/usr/bin/env Rscript

dir.create(out_dir)

mapping = read.csv(map_file, row.names=1, check.names=F)

gpr_files = list.files(gpr_dir, full.names=T)
for (k in 1:length(gpr_files)) {
	gpr_file = gpr_files[k]	
	simple_file_name = unlist(sapply(colnames(mapping), function(x)grep(paste0(x,'-'), gpr_file, value=T)))
	## find the header line
	head_row = grep('^\\"Block', scan(gpr_file, n=50, what='character', sep='\n'), value=F)
	gpr_data = read.table(gpr_file, header=T, sep='\t', check.names=F, as.is=T, skip=head_row-1)
	
	if (!is.null(ID_to_change)) {
		changeIndex = grep(ID_to_change, gpr_data$ID, value=F)
		newID = paste(gpr_data$ID[changeIndex], gpr_data$Name[changeIndex], sep='_')
		gpr_data$ID[changeIndex] = newID
	}

	gpr_split_data = split(gpr_data, f=gpr_data$Block)
	
	for (j in 1:length(gpr_split_data)){
		blockName = names(gpr_split_data)[j]	
		blockData = gpr_split_data[[j]]
		write.table(blockData, file=file.path(out_dir, paste0(mapping[blockName,names(simple_file_name)], '.txt')), sep='\t', quote=F, row.names=F)
	}
	
}



setwd("D:/Projet Safir/Data SAFIR/Safir Output/Safir SegmentTables")
library(annotate)

census <- read.csv("D:/Projet Safir/Data SAFIR/Arrays Infos/CensusTable_Annot_FC.txt", header = T, skip = 0, sep = "\t")

segtable <- read.table("S008_252206016454_1_3_09-09-2011_Agilent_SegmentsTable.txt", header = T, sep = "\t")
Safir.Id = "S008"
index <- which(abs(segtable$seg.mean)>=1)
segtable[index,]
table(segtable$chrom[index])

Append = FALSE
if(length(index)>1) Append = TRUE
for(i in index){
	start <- segtable$loc.start[i]
	end <- segtable$loc.end[i]
	census.index <- which(	(census$Genom.start>=start & census$Genom.end<=end) |
					(census$Genom.end>=start & census$Genom.start<=end))
	census.index <- unique(census.index)
	if(length(census.index)>0){
		genlist <- as.data.frame(census$GeneId[census.index])
		seg.value <- segtable$seg.mean[i]
		chr <- segtable$chrom[i]
		htable <- htmlpage(	genlist, filename = paste(Safir.Id, "SupplGeneList.html", sep = "_"), 
				title = paste("Chr ", chr, ". Segment loc: ", start, "-", end, ". Value (Log2R) =", seg.value, " (Cancer Gene Census)", sep = ""), 
				othernames = as.data.frame(census[census.index,1:9]), table.center = T,
				table.head = c("GeneId", colnames(census)[1:9]), repository = list("en"), append = Append,
				col = rep(c("royalblue1", "seagreen1"), length(genlist)))
		}
	}



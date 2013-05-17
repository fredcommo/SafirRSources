LoadSegTable <- function(Root = "E"){

setwd(paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir SegmentTables", sep = ""))
require(tcltk)

	fileName <- tclvalue(tkgetOpenFile()) # Very simple, isn't it?
	if (!nchar(fileName)) {
		tkmessageBox(message = "No file was selected!")
		}
	else {
		tkmessageBox(message = paste("The file selected was", fileName))
		}

	fnames <- unlist(strsplit(fileName, "/"))
	nnames <- length(fnames)

	path <- fnames[1]
	for(i in 2:(nnames-1)) path <- paste(path, fnames[i], sep = "/")
	Folder <- fnames[nnames-1]
	FileName <- fnames[nnames]
	seg.table <- read.csv(FileName, header = T, sep = "\t")
	Ids <- unlist(strsplit(FileName, "_"))
	SafirId <- Ids[1]
	BarCode <- paste(Ids[2], Ids[3], Ids[4], sep = "_")
	analysisDate <- Ids[5]
	cgh <- list(Safir.Id = SafirId, BarCode = BarCode, analysisDate = analysisDate)
	cat(FileName, "loaded from", path, "\n")
	return(list(seg.table = seg.table, cgh = cgh))
}
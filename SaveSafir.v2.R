SaveSafir.v2 <- function(Root = "E", Infos = T, Plot = TRUE, Greq = TRUE, SegTab = TRUE, FullTab = T, Tech = NA){

Main.Path = paste(Root, ":/Projet Safir/Data Safir/Safir Output/", sep = "")

pathPlot <- paste(Main.Path, "Safir Plots/", sep = "")
pathGeneOfInt <- paste(Main.Path, "Safir GenesOfInterest/", sep = "")
pathSegTab <- paste(Main.Path, "Safir SegmentTables/", sep = "")
pathFullTab <- paste(Main.Path, "Safir FullOutputs/", sep = "")
pathInfos <-  paste(Main.Path, "Safir ScriptInfos/", sep = "")

if(Infos){
	write.table(cgh.infos, paste(pathInfos, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_ScriptInfos.txt", sep = ""), sep = "\t", row.names = F)
	cat("Analysis infos saved:\n", paste(pathInfos, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_ScriptInfos", sep = ""), "\n\n")
	}
if(Plot){
	savePlot(filename = paste(pathPlot, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_GenomicPlot", sep = ""), type = "png")
	cat("Plot saved:\n", paste(pathPlot, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_GenomicPlot", sep = ""), "\n\n")
	}
if(Greq){
	write.table(greq, paste(pathGeneOfInt, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_GenesOfInt.xls", sep = ""), sep = "\t", row.names = F)
	cat("Genes of interest table saved:\n",  paste(pathGeneOfInt, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_GenesOfInt.xls", sep = ""), "\n\n")
	}
if(SegTab){
	write.table(seg.table, paste(pathSegTab, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_SegmentsTable.txt", sep = ""), sep = "\t", row.names = F)
	cat("Segmentation table saved:\n", paste(pathSegTab, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_SegmentsTable.txt", sep = ""), "\n\n")
	}
if(FullTab){
	write.table(Full.output, paste(pathFullTab, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_FullOutput.txt", sep = ""), sep = "\t", row.names = F)
	cat("Plot saved:\n",  paste(pathFullTab, Scgh$Safir.Id, "_", Scgh$BarCode, "_", Scgh$analysisDate, "_", Tech, "_FullOutput.txt", sep = ""), "\n\n")
	}
}


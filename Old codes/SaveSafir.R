SaveSafir <- function(){

pathPlot = "E:/Projet Safir/Data Safir/Safir Output/Safir Plots/"
pathGeneOfInt = "E:/Projet Safir/Data Safir/Safir Output/Safir GenesOfInterest/"
pathSegTab = "E:/Projet Safir/Data Safir/Safir Output/Safir SegmentTables/"
pathFullTab = "E:/Projet Safir/Data Safir/Safir Output/Safir FullOutputs/"

savePlot(filename = paste(pathPlot, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_GenomicPlot", sep = ""), type = "png")
cat("Plot saved:\n", paste(pathPlot, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_GenomicPlot", sep = ""), "\n\n")

write.table(greq, paste(pathGeneOfInt, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_GenesOfInt.xls", sep = ""), sep = "\t", row.names = F)
cat("Genes of interest table saved:\n",  paste(pathGeneOfInt, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_GenesOfInt.xls", sep = ""), "\n\n")

write.table(seg.table, paste(pathSegTab, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_Affy_SegmentsTable.txt", sep = ""), sep = "\t")
cat("Segmentation table saved:\n", paste(pathSegTab, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_Affy_SegmentsTable.txt", sep = ""), "\n\n")

write.table(Full.output, paste(pathFullTab, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_Affy_FullOutput.txt", sep = ""), sep = "\t")
cat("Plot saved:\n",  paste(pathFullTab, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_Affy_FullOutput.txt", sep = ""), "\n\n")
}


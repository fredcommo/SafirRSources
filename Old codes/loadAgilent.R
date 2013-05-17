
loadAgilent <- function(Folder, FileName, Safir.Id){

require(limma)
require(DNAcopy)
require(preprocessCore)
require(AMORE)
require(affy)
require(robustbase)
require(mclust)

source("E:\\Stats\\Doc R\\Scripts R\\SAFIR\\GeneRequest.R")
op <- par(no.readonly = T)

plot.new()

path <- "E:/Projet Safir/Data SAFIR/Safir Input"
setwd(paste(path, Folder, "",sep = "/"))

##########################################################################################

	# Step1: Chargement des données
cat("Step1.a: Collecte des informations mArray:")
a.info <- try(read.csv(paste(FileName, ".txt", sep = ""), header = F, fill = T, skip = 1, nrows = 8, sep = "\t"), silent = T)

{
if(class(a.info)!="try-error") cat("\nCGH data have been correctly downloaded:", "\n\tFolder: ", getwd(), "\n\tFileName:", FileName, "\n\tSafir.Id:", Safir.Id, "\n")
else stop("\t",FileName, ": No such file on directory (glup!). Please fix it!\n")
}

ArraySet = "Agilent-022060 SurePrint G3 Human CGH Microarray 4x180K"
Protocol = as.character(a.info[2,2])
ScanName = as.character(a.info[2,4])
ScanDate = as.character(a.info[2,6])
ScanDate <- strsplit(ScanDate, " ")[[1]][1]
GridName = as.character(a.info[2,11])
BarCode = as.character(a.info[2,27])
BarCode <- strsplit(BarCode, "_")[[1]][2]
LabId = as.character(a.info[2,31])
FEversion = as.character(a.info[2,33])
analysisDate <- format(Sys.Date(), "%m-%d-%Y")
 
cat("\nStep1.b: Chargement des données mArray:")
cgh <- read.csv(paste(FileName, ".txt", sep = ""), header = T, skip = 9, sep = "\t")
cgh <- cgh[which(substr(cgh$ProbeName, 1, 1)== "A"),c(7:8, 18:19, 24:27, 39:40)]
cat("\nCGH data have been correctly downloaded\n")

		# table d'annotation Agilent (+ positions génomiques hg19)
cat("\nStep1.c: Chargement des tables d'annotation complémentaires:")
agil.hg19 <- try(read.table("E:/Projet Safir/Data SAFIR/Arrays Infos/022060_4x180K_hg19_20110628_GenomePos&GC_FC.txt", header = T, sep = "\t"), silent = T)
{
if(class(agil.hg19)!="try-error") cat("\nAgilent-hg19 information file have been correctly downloaded\n")
else stop("\tNo Agilent hg19 on this directory (glup!). Please fix it!\n")
}

agil.hg19 <- agil.hg19[order(agil.hg19$ProbeName),]
# hg19.info <- try(read.table("E:/Projet Safir/Data SAFIR/Arrays Infos/human.chrom.info.hg19.FC.txt", header = T, sep = "\t"))
# if(class(hg19.info)=="try-error") stop("\thg19.info: No such file on directory (glup!). Please fix it!\n")


##########################################################################################

	# Step2: Identification & suppression des flags
cat("\nStep2: suppression des flags\n")
flags.info <- cgh[,5:10]
flags <- which(	flags.info[,1]== 1 | 		# 1 = gTsSaturated invalide
		flags.info[,2]== 1 | 			# 1 = rTsSaturated invalide
		flags.info[,3]== 1 | 			# 1 = gIsFeatureNonUnifOL invalide
		flags.info[,4]== 1 | 			# 1 = rIsFeatureNonUnifOL invalide
		flags.info[,5]== 0 | 			# 0 = gIsWellAboveBG invalide
		flags.info[,6]== 0)			# 0 = rIsWellAboveBG invalide

cgh <- cgh[-flags,]

##########################################################################################

	# Step3: suppression des sondes dupliquées
cat("\nStep3: Suppression des sondes dupliquées\n")
hg19.duplic <- which(duplicated(agil.hg19$ProbeName))
if(length(hg19.duplic>0)) agil.hg19 <- agil.hg19[-hg19.duplic,]
rownames(agil.hg19) <- seq(1, nrow(agil.hg19))

cgh <- cgh[order(cgh$ProbeName),]
cgh.duplic <- which(duplicated(cgh$ProbeName))
if(length(cgh.duplic>0)) cgh <- cgh[-cgh.duplic,]
rownames(cgh) <- seq(1, nrow(cgh))

##########################################################################################

	# Step4: réannot selon hg19
cat("\nStep4: Calcul des positions génomiques selon hg19\n")

hg19 <- agil.hg19[which(is.element(agil.hg19$ProbeName, cgh$ProbeName)), ]
cgh.hg19 <- cgh[which(is.element(cgh$ProbeName, agil.hg19$ProbeName)), ]

verif <- ifelse(as.character(cgh.hg19$ProbeName) == as.character(hg19$ProbeName), "ok", "error")
n.error <- length(which(verif == "error"))

{
if(n.error == 0) {
	cgh <- cbind.data.frame(cgh.hg19[,1:2], hg19, cgh.hg19[,-c(1:2)])
	cat("\n\nHow lucky you are :-)\n\n\n")
	}
else stop("\tStep4: Errors on hg19 annotations")
}
colnames(cgh)[4:6] <- c("ChrNum", "ChrStart", "ChrEnd")
cgh <- cgh[order(cgh$ChrNum, cgh$ChrStart),]

cat("\nSafir.Id:", Safir.Id, "\nBarCode:", BarCode, "\nArraySet:", ArraySet, "\nProtocol:", Protocol, "\nScanName:", ScanName, "\nScanDate:", ScanDate, "\nLabId:", LabId, "\nFEversion:", FEversion, "\nAnalysis Date:", analysisDate, "\n")
return(list(cgh = cgh, Safir.Id = Safir.Id, BarCode = BarCode, ArraySet = ArraySet, Protocol = Protocol, ScanName = ScanName, ScanDate = ScanDate, LabId = LabId, FEversion = FEversion, analysisDate = analysisDate))
}
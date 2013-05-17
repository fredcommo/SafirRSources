# CGH Robject TestSuite BlackBox

Root = "/Users/fredcommo/Documents/Projet Safir/"


# Test Affy workflow
setwd(paste(Root, 'Safir R Sources/', sep ='')); source('SourceCode.R'); source("loadCGHObj.01.R")

f = getFile(); f
a = getAnnot(f); a
info = getAffyInfo(f); info
X = getAffyData(f); head(X)
dim(X)
X2 <- nonDuplic(X); head(X2)
dim(X2)
X3 = addAnnot(X2, a$platform); head(X3)
dim(X3)


# Test Agilent workflow
setwd(paste(Root, 'Safir R Sources/', sep ='')); source('SourceCode.R'); source("loadCGHObj.01.R")

f = getFile(); f
a = getAnnot(f); a
info = getAgilentInfo(f); info
X = getAgilentData(f); head(X)
flags = flagList(X)
length(flags)
X$gMedianSignal[flags] = X$rMedianSignal[flags] = NA
head(X); dim(X)
table(X$ChrNum)

keepCol = c('ProbeName', 'ChrNum', 'ChrStart', 'gMedianSignal', 'rMedianSignal')
X <- X[ ,colnames(X) %in% keepCol]; head(X); dim(X)
system.time(X2 <- nonDuplic(X)); head(X2); dim(X2)
X3 = addAnnot(X2, a$platform); head(X3); dim(X3)


# Test global
setwd(paste(Root, 'Safir R Sources/', sep ='')); source('SourceCode.R'); source("loadCGHObj.01.R")

myObj = loadCGHObj.01()
getInfo(myObj)
getInfo(myObj, 'sampleId')
getInfo(myObj, 'fileName')

X = getCNset(myObj)
head(X)

setwd(paste(Root, 'Safir R Sources/', sep ='')); source('SourceCode.R'); source("NormCGHObj.R")
X2 = CyAdjust(X)
X3 = GCadjust(X2)


systName = as.character(X$SystematicName)
testAnnot = substr(systName, 1, 3)
notChr = which(testAnnot != 'chr')
systName[notChr] = 'chrNA:NA-NA'
splitNames = unlist(strsplit(systName, ":"))
Chr = splitNames[seq(1, length(splitNames), by = 2)]
ChrNum = substr(Chr,4,10)
randomProbes = which(nchar(ChrNum)>2)
ChrNum[randomProbes] = NA
ChrNum[which(ChrNum=='M')] = NA
ChrNum[which(ChrNum=='X')] = 23
ChrNum[which(ChrNum=='Y')] = 24
ChrNum = as.numeric(ChrNum)

positions = splitNames[seq(2, length(splitNames), by = 2)]
positions[randomProbes] = 'NA-NA'
splitPos = unlist(strsplit(positions, '-'))
start = as.numeric(splitPos[seq(1, length(splitPos), by = 2)])
end = as.numeric(splitPos[seq(2, length(splitPos), by = 2)])


dim(X); length(ChrNum); length(start); length(end)

Chr = cbind.data.frame(splitNames[seq(1, length(splitNames), by = 2)], splitNames[seq(2, length(splitNames), by = 2)])
Chr[1:10]
table(Chr)

testChr = substr(Chr[,1], 1, 3)
err = which(testChr != 'chr')
Chr[err[1:5],]

# Duplic function

nonDuplic <- function(X){
	
	'
	Suppress the duplicated probes in a given array matrix, and return the curated matrix.
	The column containing the probeIds have to be named $ProbeName
	'

	duplicIndex <- which(duplicated(X$ProbeName))
	if(length(duplicIndex>0)) {
		duplicProbes <- as.character(unique(X$ProbeName[duplicIndex]))
		sub = subset(X, X$ProbeName %in% duplicProbes, select = c(ProbeName, gMedianSignal, rMedianSignal))
		for(probe in duplicProbes){
			tmp = subset(sub, sub$ProbeName == probe, select = c(gMedianSignal, rMedianSignal))
			m = sapply(tmp, median, na.rm = TRUE)
			X$gMedianSignal[X$ProbeName == probe] = m[1]
			X$rMedianSignal[X$ProbeName == probe] = m[2]
			}
		X <- X[-duplicIndex,]
		}
	rownames(X) <- seq(1, nrow(X))
	return(X)
}


# Annot snp6_GCpercent with genomic positions

Path = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos'
gcFile = '/022060_4x180K_hg19_20110628_GCpercent.txt'
gcPercent <- read.csv(paste(Path, gcFile, sep = ''), header = TRUE, sep = '\t')
gcPercent = gcPercent[order(gcPercent$ChrNum, gcPercent$ChrStart),]
head(gcPercent)

chrLen <- read.table('/Users/fredcommo/Documents/Projet Safir/Arrays Infos/human.chrom.info.hg19.FC.txt', header = TRUE, sep = "\t")
cumLen <- cumsum(as.numeric(chrLen$length))
cumLen <- c(0, cumLen[-length(cumLen)])

gPosition = c()
for (i in 1:24)
	gPosition <- c(gPosition, rep(cumLen[i], length(which(gcPercent$ChrNum==i))))

if(any(is.na(gcPercent$ChrNum))){
	nNA = length(which(is.na(gcPercent$ChrNum)))
	gPosition <- c(gPosition, rep(NA, nNA))
	}
	
gPosition <- gcPercent$ChrStart + gPosition

end <- c()
for (i in 1:24)
	end <- c(end, max(gPosition[which(gcPercent$ChrNum==i)]))

plot(gcPercent$ChrNum, gPosition, cex = 0.1)
abline(h = end, col = 'red', lwd = 0.1)

overlap <- c()
chrList <- unique(gcPercent$ChrNum)
for(i in chrList)
	if(all(c(i-1, i) %in% chrList)){
		last = max(gPosition[which(gcPercent$ChrNum == i-1)], na.rm = TRUE)
		first = min(gPosition[which(gcPercent$ChrNum == i)], na.rm = TRUE)
		cat(i-1, 'Vs.', i, '\tLast:', last, '\tFirst:', first, '\n')
		overlap <- c(overlap, last>first)
		}
overlap

gcPercent = cbind.data.frame(gcPercent[,c(1, 3:5)], genomicPos = gPosition, GCpercent = gcPercent[,7])
head(gcPercent)

write.table(gcPercent, paste(Path, gcFile, sep = ''), row.names = F, sep = '\t')

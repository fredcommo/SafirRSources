EditInfos <- function(cgh){

	n.infos <- length(cgh)
	cgh.infos <- c()

	for(i in 2:n.infos){
		tmp <- paste(names(cgh[i]), as.character(cgh[[i]]), sep = "=")
		cgh.infos <- c(cgh.infos, tmp)
		}
	cgh.infos <- as.data.frame(cgh.infos)
	return(cgh.infos)
}
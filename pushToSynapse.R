#######################
# Push FE files in a synapse directory
#######################

# Install synapseRclient
source('http://depot.sagebase.org/CRAN.R')
pkgInstall("synapseClient")

require(synapseClient)
require(foreach)
require(iterators)

# Use your login and password here.
synapseLogin('myName@email.org', 'myPwd')
# Indicate the path to your directory here
setwd("/myPath/toFiles")

listFiles = list.files()
foreach (file = iter(listFiles)) %do% {
  newData <- Data(list(name = file, parentId = 'syn1715847')) 	# create a child Id
  newData <- addFile(newData, file)											# add the file
  newData <- storeEntity(newData)											# push in synapse
}


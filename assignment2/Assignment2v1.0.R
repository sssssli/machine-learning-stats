#Q1A
#setwd("/Users/siyangli/Documents/Grad1/Statistical modeling and R/")
#the above is where I saved "ass2-simpleclust.txt"

simpleData <- read.table ("ass2-simpleclust.txt", header=FALSE, sep = "",quote="")
simpleData
#The data are x and y coordinates for one single point, stored as a dataframe
names (simpleData) <- c("x", "y")
getXcoords <- function (df){ #input is a data.frame with numbers 
	return (as.vector(unlist(df$x)))
}

getYcoords <- function (df){
	return (as.vector(unlist(df$y)))
}

simpleXcoords <- getXcoords (simpleData)
simpleYcoords <- getYcoords (simpleData)
plot (simpleXcoords, simpleYcoords, main ="Data from 'ass2-simpleclust.txt'" )


#Q1B
q1distMat <- dist (simpleData, method = "euclidean", diag=FALSE, upper= FALSE, p=2) #p=power of Minkowski distance, p=2 when it is Euclidean
q1completeClust <- hclust (q1distMat) #the default agglomeration method used is "complete" and that the given distance matrix is based on distances between single points (not clusters)
par () #to visualise current graphical parameters
par (cex.main="1.1") #changed from 
plot (q1completeClust,main = "Cluster dendrogram of the distance matrix from 'ass2_simpleclust.txt' by \n'complete' agglomeration method", xlab = "Distance matrix of 'ass2_simpleclust.txt'")

#Q1C
q1kmeans <- kmeans (simpleData, 2, iter.max=10, nstart =1, algorithm = "Lloyd")
plot (simpleXcoords, simpleYcoords, main ="Data from 'ass2-simpleclust.txt' with \n centres obtained from k-means" ) #need to plot the graph again before adding points onto the plot
points (q1kmeans$centers, pch=16, col=12)
q1kmeansLegend <- c("ass2-simpleclust.txt data", "k-means centres")

legend ("bottomright", q1kmeansLegend, col=c(par("col"), 12), pch=c(par("pch"), 16), title="Legend")

#Q1D
#Install mclust package under CRAN binary

library (mclust)
q1Mclust <- Mclust (simpleData)
q1Mclust
plot (q1Mclust, data=simpleData)

 # best model: elliposidal, equal variance with 2 components
 # plot (q1Mclust)

q1MclustBIC <- mclustBIC (simpleData)
q1MclustSummary <- summary (q1MclustBIC, data=simpleData)
q1MclustSummary
 
 
#Q2A
hardData <- read.table ("ass2-hardclust.txt")
names (hardData) <- c("x", "y")
hardXcoords <- getXcoords (hardData)
hardYcoords <- getYcoords (hardData)
plot (hardXcoords, hardYcoords, main ="Data from 'ass2-hardclust.txt'" )


#Q2B
q2distMat <- dist (hardData, method = "euclidean", diag=FALSE, upper= TRUE, p=2) #p=power of Minkowski distance, p=2 when it is Euclidean. dist() calculates between ROWS 
q2completeClust <- hclust (q2distMat) 
plot (q2completeClust,main = "Cluster dendrogram of the distance matrix from 'ass2_hardclust' by \n'complete' agglomeration method")
#need to fix the dendrogram 

#heat map
q2distMatMA <- as.matrix (q2distMat)
heatmap (q2distMatMA, Rowv=as.dendrogram(q2completeClust), Colv="Rowv", symm=TRUE)


#Q2C

#K means with x number of centers 

kmeansAndPlot <- function (data, numberOfCenters, iter){ #data is in the form of a data.frame, no need to input the algorithm method because we will always be using Lloyd, no need to change the number of iterations (always 10), nstart is always 1. 
	kmeansResults <- kmeans (data, numberOfCenters, iter.max=iter, nstart=1, algorithm="Lloyd") 
	#Calls other functions
	dataXcoords <- getXcoords (data)
	dataYcoords <- getYcoords (data)
	plotCentersKmeans (dataXcoords, dataYcoords, kmeansResults$centers)
	plotCenterStdev (kmeansResults, numberOfCenters)
}

plotCentersKmeans <- function (Xcoords, Ycoords, centersCoords) {
	plot (Xcoords, Ycoords, main = paste("Data from 'ass2-hardclust.txt' \n with k-means clustering algorithm (", nrow(centersCoords), " centers)", sep="" ))
	points (centersCoords, pch=16, col=34)
	kmeansLegend <- c("data points", "k-means centres")
	#legend ("bottomright", kmeansLegend, col=c(par("col"), 34), pch=c(par("pch"), 16), title="Legend")
}


#plotting the standard deviations requires knowing the standard deviation of each cluster. The kmeans () function gives the sum of squares by cluster. The standard deviation of each cluster is sqrt (withinss for each cluster/(n-1)), where n= the number of samples in that cluster. And to plot the standard deviation, we can just plot a circle.

#download the plotrix library
library (plotrix)

plotCenterStdev <- function (kmeansResults,numberOfCenters){ #draw.circle only takes x and y coords for the centers, so we need to get this from kmeansResults 
	for (i in 1:numberOfCenters){
		circleX <- kmeansResults$centers [i,1]
		circleY <- kmeansResults$centers [i,2]
		stdev <- sqrt (kmeansResults$withinss [i]/(kmeansResults$size[i]-1))
		draw.circle (circleX, circleY, stdev, border="red", lwd=1.2)
	}
	kmeansLegend <- c("data points", "k-means centres", "standard deviation")
	legend ("bottomright", kmeansLegend, col=c(par("col"), 34, "red"), pch=c(par("pch"), 16), title="Legend")
	
}
kmeansAndPlot (hardData, 2)
kmeansAndPlot (hardData, 4)
kmeansAndPlot (hardData, 10)


#Q2D
q2Mclust <- Mclust (hardData)
q2Mclust 

# best model: elliposidal, equal variance with 2 components

plot (q2Mclust, data=hardData)

q2MclustBIC <- mclustBIC (hardData) 
q2MclustSummary <- summary (q2MclustBIC, data=hardData)
q2MclustSummary


#Q3A need to do hclust() and Mclust() and time the R commands on randomly drawn samples from our list
MicroArray <- read.table ("microarray_data.txt", header=T, sep="\t")
head (MicroArray) #Data is arranged by columns (each experiment is a new column) 
matrixMA <- as.matrix (MicroArray)
matrixMA <- matrixMA [rowSums(is.na(matrixMA))==0,] #removes all genes with missing datapoints
MicroArray <- as.data.frame (matrixMA)

#randomly drawing samples
totalGenes <- nrow (MicroArray) 
colMicroArray <- length (MicroArray)
sampling <- function (sampleNumber){
	geneIndices <- sample (totalGenes, sampleNumber)
	randomSamples <- sapply (1:sampleNumber, function(x)matrixMA[geneIndices[x],])
	randomSamples <- t(randomSamples) #transposes the matrix, since dist only performs between rows 
	randomSamplesNoName <- randomSamples [,2:colMicroArray]
	class (randomSamplesNoName) <- "numeric"
	return (randomSamplesNoName)
	
}


timeCalculation <- function (randomGenes, clustMethod){
	if (clustMethod =="hclust"){
		q3distMat <- dist (randomGenes, method = "euclidean", diag=FALSE, upper= FALSE, p=2)
		q3completeClust <- hclust(q3distMat)
		
	}else if (clustMethod =="Mclust"){
		Mclust (randomGenes)
		
	}
}
random200 <- sampling (200)
timeH200 <- system.time (timeCalculation(random200, "hclust")) #system CPU time is the second number, so timeH200 [2] would call this time 
timeM200 <- system.time (timeCalculation(random200, "Mclust"))

random500 <- sampling (500)
timeH500 <- system.time (timeCalculation(random500, "hclust")) 
timeM500 <- system.time (timeCalculation(random500, "Mclust"))

random1000 <- sampling (1000)
timeH1000 <- system.time (timeCalculation(random1000, "hclust")) 
timeM1000 <- system.time (timeCalculation(random1000, "Mclust"))

random2000 <- sampling (2000)
timeH2000 <- system.time (timeCalculation(random2000, "hclust")) 
timeM2000 <- system.time (timeCalculation(random2000, "Mclust"))

timeH <- cbind (timeH200[3], timeH500[3], timeH1000[3], timeH2000[3], timeH4880[3])
timeM <- cbind (timeM200[3], timeM500[3], timeM1000[3], timeM2000[3], timeM4880[3])

combinedData <- rbind (timeH, timeM)
rownames (combinedData) <- cbind ("hclust()", "Mclust()")
colnames (combinedData) <- cbind (200, 500, 1000, 2000)
combinedData

barplot (combinedData, main = "Comparison of hierarchical clustering (hclust ()) and \nGaussian mixture model-based clustering (Mclust())\nfor different numbers of data points", xlab="Number of random samples drawn from 'micro_array.txt'", ylab="Elapsed time (s) for algorithm to run", legend = TRUE, beside=TRUE, args.legend = list(x="topleft"), ylim=c(0,30)) #legend=TRUE will just print the names of the rows of combined, args.legend sets the legend to topleft of the graph


#Q3B, let's keep the number of datapoints at 1000
random1000for3B <- sampling (1000)
class (random1000for3B) <- "numeric"
randomMicroArrayExp <- function (numberOfExperiments){
	columns3B <- sample (2:colMicroArray, numberOfExperiments)
	return (random1000for3B[,columns3B]) #returns all the columns that were selected under the random column numbers 
}

random1000_2exps <- randomMicroArrayExp (2)
random1000_4exps <- randomMicroArrayExp (4)
#for the 8 microarray experiments, no need to sample, since we only have 8 micro array experiments in total!
random1000_8exps <- random1000for3B

timeH2exps <- system.time (timeCalculation(random1000_2exps, "hclust"))
timeM2exps <- system.time (timeCalculation(random1000_2exps, "Mclust"))

timeH4exps <- system.time (timeCalculation(random1000_4exps, "hclust"))
timeM4exps <- system.time (timeCalculation(random1000_4exps, "Mclust"))

timeH8exps <- system.time (timeCalculation(random1000_8exps, "hclust"))
timeM8exps <- system.time (timeCalculation(random1000_8exps, "Mclust"))

timeHdimensions <- cbind (timeH2exps[3], timeH4exps[3], timeH8exps[3])
timeMdimensions <- cbind (timeM2exps[3], timeM4exps[3], timeM8exps[3])

combinedDims<- rbind (timeHdimensions, timeMdimensions)
rownames (combinedDims) <- cbind ("hclust()", "Mclust()")
colnames (combinedDims) <- cbind ("2 dimensions", "4 dimensions", "8 dimensions")
combinedDims
barplot (combinedDims, main = "Comparison of hierarchical clustering (hclust()) and \nGaussian mixture model-based clustering (Mclust())\nfor different numbers of dimensions on 1000 data points", xlab="Number of dimensions (randomly drawn) for 1000 data points from 'micro_array.txt'", ylab="Elapsed time (s) for algorithm to run", legend = TRUE, beside=TRUE, args.legend = list(x="topleft"), ylim=c(0,15)) #legend=TRUE will just print the names of the rows of combined, args.legend sets the legend to topleft of the graph


#Q3C
q3 <- cbind (as.numeric(matrixMA[,2]),as.numeric(matrixMA[,3]),as.numeric(matrixMA[,4]),as.numeric(matrixMA[,5]),as.numeric(matrixMA[,6]),as.numeric(matrixMA[,7]),as.numeric(matrixMA[,8]), as.numeric(matrixMA[,9]))  #they were somehow stored as strings in a matrix before and this is the only way to convert them 

q3mclust <- Mclust (q3)
q3mclust

# best model: ellipsoidal, equal shape with 6 components. But 6 clusters may not always be right, since the question asks for us to use kmeans (), the below tries a few different centres for kmeans ()

q3_4kmeans <- kmeans (q3, 4)
q3_5kmeans <- kmeans (q3, 5)
q3_6kmeans <- kmeans (q3, 6)
q3_7kmeans <- kmeans (q3, 7)
q3_8kmeans <- kmeans (q3, 8, iter.max=20) #no convergence in 10 iterations, changed to 20

q3_4kmeans$size
#[1] 1763  656 2453    8

q3_5kmeans$size
#[1] 1128  855 1635   38 1224

q3_6kmeans$size
#[1] 1110 1137  133 1515    8  977

q3_7kmeans$size
#[1] 1413  827  664 1010    8  839  119

q3_8kmeans$size
#[1] 1019  646  107 1072  336  666 1026    8

q3_6kmeans$centers
        [,1]       [,2]        [,3]       [,4]        [,5]       [,6]        [,7]        [,8]
1  0.29740541  0.1031982 -0.03678378 -0.1094685  0.33325225 -0.4685135  0.54597297  0.14727928
2 -0.28160070 -0.1111170 -0.13931398 -0.3055409 -0.13277045 -0.9126737 -0.07798593 -0.02496042
3  0.51187970  0.1681203  0.13563910  2.8430075  0.65526316  0.2143609  0.67857143  0.16390977
4  0.07376898 -0.0310033 -0.01849505  0.4189175  0.07679868 -0.2735776  0.21266007 -0.07854125
5  3.95500000  3.9825000  4.14250000  4.5850000  4.52750000  2.8587500  3.38500000  2.55500000
6  0.47750256  0.1504811  0.12264074  0.7052815  0.47326510  0.3149130  0.87508700  0.10855681
q3_4kmeans$centers
        [,1]        [,2]         [,3]       [,4]        [,5]       [,6]         [,7]        [,8]
1 -0.1490868 -0.06747589 -0.123028928 -0.2114804 -0.05311401 -0.8107657 -0.009699376 -0.02003971
2  0.5726372  0.17679878  0.164984756  1.4075457  0.57853659  0.3121189  0.827896341  0.12958841
3  0.2215614  0.04631064  0.007859764  0.2781981  0.24417040 -0.1677212  0.513147167  0.03843865
4  3.9550000  3.98250000  4.142500000  4.5850000  4.52750000  2.8587500  3.385000000  2.55500000
q3_5kmeans$centers
        [,1]        [,2]         [,3]       [,4]        [,5]       [,6]       [,7]        [,8]
1  0.3313209  0.13852837 -0.050735816 -0.0656117  0.37306738 -0.4205496  0.5931472  0.15371454
2  0.5114035  0.16008187  0.130327485  1.0271345  0.52759064  0.3071813  0.8664444  0.10796491
3  0.0804526 -0.03798777  0.004214067  0.4057859  0.08264832 -0.2142813  0.2579572 -0.06574312
4  1.5800000  1.22815789  1.359210526  4.3073684  1.68763158  1.1452632  1.5426316  0.84552632
5 -0.2570425 -0.11046569 -0.138063725 -0.2955229 -0.12264706 -0.8942892 -0.0640768 -0.02118464
q3_6kmeans$centers
         [,1]       [,2]        [,3]       [,4]        [,5]       [,6]        [,7]        [,8]
1  0.29740541  0.1031982 -0.03678378 -0.1094685  0.33325225 -0.4685135  0.54597297  0.14727928
2 -0.28160070 -0.1111170 -0.13931398 -0.3055409 -0.13277045 -0.9126737 -0.07798593 -0.02496042
3  0.51187970  0.1681203  0.13563910  2.8430075  0.65526316  0.2143609  0.67857143  0.16390977
4  0.07376898 -0.0310033 -0.01849505  0.4189175  0.07679868 -0.2735776  0.21266007 -0.07854125
5  3.95500000  3.9825000  4.14250000  4.5850000  4.52750000  2.8587500  3.38500000  2.55500000
6  0.47750256  0.1504811  0.12264074  0.7052815  0.47326510  0.3149130  0.87508700  0.10855681




#Cluster of 8 consistantly appears with centres >= 6. Because of the information that we were given by the question, that there is a small set of 
MAwithClust <- cbind (matrixMA, q3_6kmeans$cluster)
q3cluster <- MAwithClust [which(MAwithClust[,10]==5),]
q3cluster

q3clusterNames <- q3cluster[,1]
q3clusterNames <- as.data.frame (q3clusterNames) #needs a data frame for getGOinList

# Download org.Sc.sgd.db package
library ("org.Sc.sgd.db")
ygenes <- org.Sc.sgdGO

#Some functions from previous assignment

#Function that gets all the GO terms in a gene list

getGOinList <- function (geneList){ #geneList should be a dataframe and the first column should be entrez genes
	sampleGOIDs <- list () #empties list from previous usage
	sampleGOIDs <- sapply (1:nrow(geneList), function(x)names(ygenes[[as.character(geneList[x,1])]]))
	#geneList[i,1] looks up the entrez ID through the entire list and the names function retrieves all the GO IDs that is associated with each entrez gene. I used as.character here to convert the entrez ID into a string, as I said previously that looking up with the index number would be incorrect. Instead of doing this in a loop, I chose to use the sapply function because this is much faster (it applies the function simultaneously to all vectors in the list). 
geneList$GOIDs <- sampleGOIDs
return (geneList)
}

clusterGO <- getGOinList (q3clusterNames)
justGOsample <- unique (unlist (clusterGO$GOIDs))

#All the ORFs that are linked to a GO, this is our universe...
orfsWithGO <- mappedkeys (ygenes)
numberOfGenesUniverse <- nrow(orfsWithGO)
orfsWithGO <- as.data.frame (orfsWithGO)
universeGO <- getGOinList (orfsWithGO)
justGOuniverse <- unique (unlist(universeGO$GOIDs))


hyperTest <- function (sampleGOcounts, universeGOcounts){ 
		p <- phyper (sampleGOcounts-1, universeGOcounts, nrow(orfsWithGO), nrow(clusterGO) , lower.tail = FALSE, log.p=FALSE) #totalGO is the total number of GO terms, because we already tested the GO terms that are not in our sample list, but we still need to correct for it
		return (p)
}

checkingEveryGOwithHyper <- function (GOIDs, sampleGenes){ 
#no need for universe, always the same. The sampleGenes should have GOIDs as a header and should be a data.frame. 
	
sampleListTOI <- vector () #empties vector from previous usage
universeListTOI <- vector ()
hyperTestResults <- vector ()	
for (i in 1:length(GOIDs)){ #for every GO (vector of strings)
sampleListCounts <- sum(sapply (sampleGenes$GOIDs, function(x)GOIDs[i]%in%x)) #Use sapply to look for the GO ID in our sample list's column "GOIDs". Sapply returns boolean variables, and the sum just takes the number of terms that are TRUE in sapply. So this just returns the number of genes in the list with the GOID that we're interested in. Since we are using a boolean function, it would ignore GO terms appearing more than once under the same gene (avoiding annotation mistakes)
universeListCounts <- sum(sapply (universeGO$GOIDs, function(x)GOIDs[i]%in%x))
#The universe is always the same (to simplify things)
sampleListTOI <- append (sampleListTOI, sampleListCounts)
universeListTOI <- append (universeListTOI, universeListCounts)
hyperTestResults <- append (hyperTestResults, hyperTest(sampleListCounts, universeListCounts))
#Calls the function with these parameters. 
	
}

sampleGenesAllPvalues <- list (GOIDs=GOIDs, sample_counts=sampleListTOI, universe_counts=universeListTOI, hypergeom_pvalues=hyperTestResults) #this puts all info together into one var as R functions can't return more than 1 var
sampleGenesAllPvalues <- as.data.frame (sampleGenesAllPvalues) #converts to df, required format for following functions
return (sampleGenesAllPvalues)

}


#Significance test
significanceTest <- function (alpha, geneList){ #geneList should be in the format of the output of checkingEveryGOwithHyper
	significantIndices <- which(geneList$hypergeom_pvalues < alpha)
	significantTerms <- list(GOIDs = geneList$GOIDs[significantIndices], sample_counts = geneList$sample_counts[significantIndices], universe_counts = geneList$universe_counts[significantIndices], hypergeom_pvalues=geneList$hypergeom_pvalues[significantIndices])
	significantTerms <- as.data.frame (significantTerms)
	return (significantTerms)
}

retrieveGOterm <- function (significantTerms){ #significantList is the combined matrix from printListWithoutGOterm
	if (nrow (significantTerms)==0){
		return ("No significant terms")
	}else{
	GOterms <- sapply (1:nrow(significantTerms), function(x)as.character(Term(as.character(significantTerms[[x,1]]))))
	significantTerms <- cbind (GOterms=GOterms, significantTerms)
	significantTerms <- as.data.frame (significantTerms)
	return(significantTerms)
}
}

GOclustAllP <- checkingEveryGOwithHyper (justGOsample, clusterGO)
print (data.frame (GOclustAllP))
alpha <- 0.05
GOclustSigP <- significanceTest (alpha, GOclustAllP)
final <- retrieveGOterm (GOclustSigP)
final

#Multiple testing correction with FDR
FDRAllP <- GOclustAllP
FDRAllP$hypergeom_pvalues <-  p.adjust (FDRAllP$hypergeom_pvalues, method = "fdr")
FDRsigTerms <- significanceTest (alpha, FDRAllP)
FDRsigGOterms <- retrieveGOterm (FDRsigTerms)
names (FDRsigGOterms) <- c("GOterms", "GOIDs", "sample_counts", "universe_counts", "corrected_pvalues")
FDRsigGOterms 


#Multiple testing correction with Bonferroni
BonfAllP <- GOclustAllP
BonfAllP$hypergeom_pvalues <-  p.adjust (BonfAllP$hypergeom_pvalues, method = "bonferroni", length(justGOuniverse))
BonfsigTerms <- significanceTest (alpha, BonfAllP)
BonfsigGOterms <- retrieveGOterm (BonfsigTerms)
BonfsigGOterms 
justGOuniverse



#retrieves GO terms associated with GOIDs
#completeGOterms <- as.list(GOTERM)
retrieveGOterm <- function (significantTerms){ #significantList is the combined matrix from printListWithoutGOterm
	if (nrow (significantTerms)==0){
		return ("No significant terms")
	}else{
	GOterms <- sapply (1:nrow(significantTerms), function(x)as.character(Term(as.character(significantTerms[[x,1]]))))
	significantTerms <- cbind (GOterms=GOterms, significantTerms)
	significantTerms <- as.data.frame (significantTerms)
	return(significantTerms)
}
}








pairs (q3) #To visualize 2 dimensions of the data at the same time
#Plan: trial and error with Mclust

q3kmeans <- kmeans (q3, 2)
q3kmeans #
which (q3kmeans[[1]]==1)
# xkmeans [[1]]-> cluster assignment





q3distMat <- dist (q3, method = "euclidean", diag=FALSE, upper= FALSE, p=2)
q3completeClust <- hclust (q3distMat) #the default agglomeration method used is "complete" and that the given distance matrix is based on distances between single points (not clusters)
q3completeClust









 #######
 m1 <- kmeans(ass2hardclust, 2)
 m1[[1]]
 plot(ass2hardclust[ which(m1[[1]] == 1), ], xlim=c(-2,3), ylim=c(-2,3))
> points(ass2hardclust[ which(m1[[1]] == 2), ], col='red')
 
 > plot(ass2hardclust[ which(m2[[9]] == 1), ], xlim=c(-2,3), ylim=c(-2,3))
> points(ass2hardclust[ which(m2[[9]] == 2), ], col='red')



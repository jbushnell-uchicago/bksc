rm(list=ls())

distFunc <- function(aX, aY, listOfXs, listOfYs) {
  distances <- sqrt( (aX - listOfXs)^2 + (aY - listOfYs)^2 ) 
  return (distances)
}


balancedKMeansWithPredefinedSubclusters <- function(coordsX, coordsY, subclusterAssignments, k = 2, subclusterCostsFilepath = '~/subclusterCosts.txt', clusterDistancesFilepath = '~/clusterDists.txt') {
  #Find the k-means cluster centers
  clustersNeeded <- k
  clusterCenters <- kmeans(cbind(coordsX, coordsY), clustersNeeded)
  clusterCenters <- clusterCenters$centers
  
  #Find the centers of the subclusters
  subclusterCentersX <- aggregate(coordsX~subclusterAssignments, FUN=mean)
  subclusterCentersY <- aggregate(coordsY~subclusterAssignments, FUN=mean)
  
  #Find order of nearness of subclusters to each cluster center
  nearnessOfSubclustersToClusters <- c()
  for (rowIndex in 1:nrow(clusterCenters)) {
    aX <- clusterCenters[rowIndex,1]
    aY <- clusterCenters[rowIndex,2]
    aDistances <- distFunc(aX, aY, subclusterCentersX$coordsX,subclusterCentersY$coordsY)
    orderResult <- order(aDistances)
    nearestSubclusterOrder <- rep(0, length(orderResult))
    for (orderIndex in 1:length(orderResult)) {
      aOrder <- orderResult[orderIndex]
      nearestSubclusterOrder[aOrder] <- orderIndex
    }
    nearnessOfSubclustersToClusters <- rbind(nearnessOfSubclustersToClusters, nearestSubclusterOrder)
  }
  
  #Now assign that order of nearness to each of point 
  #based on their nearness to cluster centers and membership of the ASM
  subClusterCosts <- c()
  for (rowIndex in 1:nrow(nearnessOfSubclustersToClusters)) {
    costsForCluster <- rep(0, length(subclusterAssignments))
    for (anIndex in 1:length(subclusterAssignments)) {
      aSubclusterAssignment <- subclusterAssignments[anIndex]
      costsForCluster[anIndex] <- nearnessOfSubclustersToClusters[rowIndex, aSubclusterAssignment]  
    }
    subClusterCosts <- rbind(subClusterCosts, costsForCluster)
  }
  
  #Now calculate each points distance to a cluster center
  clusterDists <- c()
  for (rowIndex in 1:nrow(subClusterCosts)) {
    distancesForACluster <- distFunc(clusterCenters[rowIndex, 1], clusterCenters[rowIndex, 2], coordsX, coordsY)
    clusterDists <- rbind(clusterDists, distancesForACluster)
  }
  
  #WRITE DATA TO USE IN JULIA OPTIMIZATION
  #subclusterCosts is the subclusterCosts variable in Julia
  #clusterDists is the dists variable in Julia
#  aSubclusterCostsFormatted <- cbind(rep('[', nrow(subClusterCosts)), subClusterCosts, rep(']', nrow(subClusterCosts)), c(rep(';', nrow(subClusterCosts) - 1), ''))
#  write.table(aSubclusterCostsFormatted, subclusterCostsFilepath, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
  aSubclusterCostsFormatted <- cbind( subClusterCosts )
  write.table(aSubclusterCostsFormatted, subclusterCostsFilepath, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
  
#  aDistsFormatted <- cbind(rep('[', nrow(clusterDists)), clusterDists, rep(']', nrow(clusterDists)), c(rep(';', nrow(clusterDists) - 1), ''))
#  write.table(aDistsFormatted, clusterDistancesFilepath, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
  aDistsFormatted <- cbind( clusterDists )
  write.table(aDistsFormatted, clusterDistancesFilepath, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
  
}

set.seed(123)
points1x <- rnorm(100, 3, 1)
points1y <- rnorm(100, 3, 1)
points2x <- rnorm(100, 5, 1)
points2y <- rnorm(100, 5, 1)
examplePointsX <- c(points1x, points2x)
examplePointsY <- c(points1y, points2y)
exampleSubclusters <- rep(1, 200)
exampleSubclusters[examplePointsX<3.5] <- 2
exampleSubclusters[examplePointsX>5 & examplePointsY>5] <- 5
exampleSubclusters[examplePointsX>5 & examplePointsY<5] <- 4
exampleSubclusters[examplePointsX<=5 & examplePointsY<=3] <- 3
exampleClustersNeeded <- 2



balancedKMeansWithPredefinedSubclusters(coordsX = examplePointsX, coordsY = examplePointsY,
                                        subclusterAssignments = exampleSubclusters,
                                        k = exampleClustersNeeded,
                                        subclusterCostsFilepath = '~/subclusterCosts.txt',
                                        clusterDistancesFilepath = '~/clusterDists.txt')
#now process in Julia
#...

#read in results from Julia
assignedClusters <- read.csv('/home/james/test.txt',sep=",",header=FALSE)
plot(examplePointsX,examplePointsY,col=apply(X = assignedClusters, MARGIN = 2, FUN = which.max), pch = exampleSubclusters)

table(apply(X = assignedClusters, MARGIN = 2, FUN = which.max), exampleSubclusters)
table(apply(X = assignedClusters, MARGIN = 2, FUN = which.max))






pointsX <- HalstedOff$Latitude
pointsY <- HalstedOff$Longitude
subclusterAssigns <- factor(HalstedOff$LSB.ASM)
plot(pointsX, pointsY, col=subclusterAssigns)
clustersNeeded <- 3
balancedKMeansWithPredefinedSubclusters(coordsX = pointsX, coordsY = pointsY,
                                        subclusterAssignments = subclusterAssigns,
                                        k = clustersNeeded,
                                        subclusterCostsFilepath = '~/subclusterCosts.txt',
                                        clusterDistancesFilepath = '~/clusterDists.txt')
#now process in Julia
#...
#after Julia, read in results from Julia
assignedClusters <- read.csv('/home/james/test.txt',sep=",",header=FALSE)
plot(pointsX,pointsY,col=apply(X = assignedClusters, MARGIN = 2, FUN = which.max), pch = as.integer(factor(subclusterAssigns)))

table(apply(X = assignedClusters, MARGIN = 2, FUN = which.max), subclusterAssigns)
table(apply(X = assignedClusters, MARGIN = 2, FUN = which.max))

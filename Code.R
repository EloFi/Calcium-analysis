###This code is used to process data for the article "Spatiotemporal reorganization of corticostriatal network dynamics encodes motor skill learning " 
#Copyright (C) 2016-2020 F. Appaix, Univ. Grenoble Alpes, INSERM U1216, Grenoble Institut Neurosciences, Grenoble, France.
    
    #This program is free software: you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation, either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.

    #You should have received a copy of the GNU General Public License
    #along with this program. If not, see: https://www.gnu.org/licenses/.


###R Packages
library(zoo)
library(signal)
library(grDevices)
library(sp)
library(factoextra)
library(cluster)
library(fpc)

library(RColorBrewer)
library(calibrate)
library(deldir)
library(gplots)
library(dplyr)
library(raster)
library(memisc)
library(mgcv)

###Open dataset
Results <- read.table("C:/Users/Dell/Desktop/Results.txt", header = TRUE, sep="\t")
Resultsxy <- read.table("C:/Users/Dell/Desktop/Resultsxy.txt", header = TRUE, sep="\t")

###Time for 15.3Hz imaging
Fq <- 15.3
dt <- (1/Fq)
Time <- (Results$X - 1) * dt * 1000

####Remove X column
Results$X <- NULL
Resultsxy$X.1 <- NULL

####Window size for rolling window 
windowSizeF0 <- Fq * 10

#Total number of Stimulations
Nb = 7
WindowStim <- length(Time)/Nb

####Create a 10sec rolling window with 50% of lowest values
F0 <- as.data.frame(rollapply(Results, width=windowSizeF0, FUN = function(z) mean(z[z < quantile(z, 0.5)]), partial=TRUE, fill=NULL, by.column=TRUE))

###DF/F0
DF <- (Results - F0)/F0
DF <- as.data.frame(mapply(FUN=function(x) ifelse(!is.na(x), x, 0), DF))

###Savitzky Golay filter 
windowSizeSG <- 7
DFfilt <- as.data.frame(mapply(FUN=function(v) sgolayfilt(v, p=3, n=windowSizeSG), DF))

####Find 5 peaks of Maximum Amplitude : Y=Amplmax & X=Peaks(Index); Threshold: 2*SD, from a part of the non-filtered signal
FindPeaksRes <- as.matrix(rollapply(Results, width=WindowStim, FUN=function(x) max(x), by=WindowStim))

row <- FindPeaksRes
for (i in 1:ncol(Results))
{
  for (j in 1:nrow(FindPeaksRes)){
    z <- which(Results[ ,i] == FindPeaksRes[j,i])
    if (length(z) > 1){
      z <- z[[1]]
      row[j,i] <- z
    }else{row[j,i] <- which(Results[ ,i] == FindPeaksRes[j,i])}
  }
}

n <- data.frame(matrix(row[1, ] - (WindowStim/2)))
for (i in 2:(Nb)){
  n[,i] <- (row[i, ] - (WindowStim/2))
}
colnames(n)[1] <- 'V1' 
for (i in 1:ncol(n))
{
  for (j in 1:nrow(n)){
    n[j,i][which(n[j, i] < 0)] <- 0
  }}

WindowStimdf <- data.frame(matrix(row[2, ] - row[1, ]))
for (i in 2:(Nb-1)){
  WindowStimdf[,i] <- (row[(i+1), ] - row[i, ])
}
colnames(WindowStimdf)[1] <- 'V1'

###Stimulations kept
nStim <- 5

WindowStimdf1 <- mean(WindowStimdf[ ,1], na.rm=T)
row7 <- mean(row[[7]], na.rm=T)

if(Nb>=7){
  FindPeaks <- as.matrix(rollapply(DFfilt[(WindowStimdf1:row7), ], width=WindowStim, FUN=function(x) max(x), by=WindowStim))
  thres <- as.matrix(mapply(FUN=function (x) 2 * sd(x), DF[(WindowStimdf1:row7), ]))
}else{
  FindPeaks <- as.matrix(rollapply(DFfilt[(WindowStimdf1:length(Time)), ], width=WindowStim, FUN=function(x) max(x), by=WindowStim))
  thres <- as.matrix(mapply(FUN=function (x) 2 * sd(x), DF[(WindowStimdf1:length(Time)), ]))
}

Peaks <- FindPeaks
for (i in 1:ncol(DFfilt))
{
  for (j in 1:nrow(FindPeaks)){
    z <- FindPeaks[j,i]
    if(z >= thres[[i]]){
      Peaks[j,i] <- as.matrix(which(DFfilt[ ,i] == z))
    }else{
      FindPeaks[j,i] <- 0
      Peaks[j,i] <- NA
    }
  }
}

for (i in 1:ncol(DFfilt))
{
  for (j in 1:(nrow(Peaks)-1)){
    Peaks[j,i][which(Peaks[j,i] < WindowStimdf1)] <- NA
    Peaks[(j+1),i][which(Peaks[(j+1),i] < (Peaks[j,i] + .5 * WindowStim))] <- NA
  }
}

for (i in 1:ncol(DFfilt))
{
  for (j in 1:nrow(FindPeaks)){
    FindPeaks[j,i][which(is.na(Peaks[j,i]))] <- NA
  }
}
Ymax <- as.data.frame(FindPeaks)
Ymax <- as.data.frame(mapply(FUN=function(x) ifelse(!is.na(x), x, 0), Ymax))

Peaksdf <- as.data.frame(Peaks)

####Find Time (Xmax) for each peak
Xmax <- data.frame()
for (i in 1:ncol(DFfilt))
{
  for (j in 1:nrow(Peaksdf)){
    if (!is.na(Peaksdf[j,i])){
      Xmax[j,i] <- lapply(FUN=function(x) Time[x], Peaksdf[j,i])
    }else{Xmax[j,i] <- 0}
  }
}
colnames(Xmax) <- colnames(DFfilt)

####Find start of response Y=Yst & X=Starts(Index)
DFfiltmean <- data.frame()
for (i in 1:ncol(DFfilt)){
  for (j in 1:Nb){
    DFfiltmean[j,i] <- as.data.frame(mean(DFfilt[(n[i,j]:(n[i,j]+(windowSizeF0/3))), i]))
  }
}
colnames(DFfiltmean) <- colnames(DFfilt)

Starts <- Peaksdf
for (i in 1:ncol(DFfilt))
{
  for (k in 1:nrow(Peaksdf)){
    if((Peaksdf[k,i] >= mean(Peaks[1, ], na.rm=T)/2) & (!is.na(Peaksdf[k,i]))){
      z <- Peaksdf[k,i]
      Startstemp <- DFfilt[z,i]
      while (Startstemp > DFfiltmean[(k+1), i])
      {
        z <- z - 1
        Startstemp <- DFfilt[z,i]
      }
      val <- DFfilt[(z-1),i]
      while (val <= Startstemp)
      {
        z <- z - 1
        Startstemp <- DFfilt[z,i]
        val <- DFfilt[(z-1),i]
      }
      Starts[k,i] <- z
    }else{Starts[k,i] <- 1}
  }}

Startsdf <- as.data.frame(Starts)

Yst <- Startsdf
for (i in 1:ncol(Startsdf)){
  for (j in 1:nrow(Startsdf)){
    if (Startsdf[j,i]!= 1){
      Yst[j,i] = DFfilt[(Startsdf[j,i]),i]
    }else{Yst[j,i] <- 0}
  }}
names(Yst) <- names(Startsdf)

####Find Time (Xst) for each start  
Xst <- as.data.frame(lapply(FUN=function(x) Time[x], Startsdf))

###Mean Amplitude Peak/cell
Ydelta <- as.data.frame(mapply('-', Ymax, Yst))
Ydelta[Ydelta == 0] <- NA
Mean <- as.matrix(apply(Ydelta, 2, mean, na.rm=TRUE))

Min <- mapply(FUN=function (x) min(x), DF)
Min <- min(Min)
Max <- mapply(FUN=function (x) max(x), DF)
Max <- max(Max)
if(Nb>=7){
  th <- as.matrix(mapply(FUN=function (x) (mean(x) + 2 * sd(x)), DF[c(WindowStim:((nStim+1) * WindowStim)), ]))#keep the same window as above, Findpeaks part.
}else{
  th <- as.matrix(mapply(FUN=function (x) (mean(x) + 2 * sd(x)), DF[c(WindowStim:(Nb * WindowStim)), ]))#keep the same window as above, Findpeaks part.
}

###Label non-MSNs
nMSN <- matrix(c(20, 63, 64, 67, 83, 85))#Put 0, if none

####Create 'Table' with (X, Y), AmplPeak (Mean), and Cells (n°)
MeanAmpl <- as.data.frame(Mean, row.names = NULL)
Table <- MeanAmpl
colnames(Table)[1] <-'MeanAmpl'
Table$X <- Resultsxy$X
Table$Y <- Resultsxy$Y
Table$Cells <- c(1:nrow(Mean))

####Add column 'Cell types' in Table
Table$Celltype <- paste0("MSN", Table[c(1:nrow(Mean)), 5])
Table[nMSN, 5] <- "nonMSN"
nonMSN <- Table[nMSN,]
nonMSN <- data.frame(nonMSN[ , 2:4])

####Add column 'Active Cells' in Table
th <- round(th, 1)
Meanr <- round(Mean, 1)
pm <- as.matrix(Meanr > th)
Table$ActCell <- as.numeric(pm)

#### % of Active cells 
ActCell <- which(Table$ActCell == 1)
NonActCell <- which(Table$ActCell == 0)
Table[NonActCell, 1] <- 0

RActCell <- (length(ActCell)/length(Resultsxy$X)) * 100 # % act cells

###Order Table by Ampl (lowest to highest)
Tableorder <- Table[order(Table[ , 1]), ]

###Create color palette with amplitude
Tableorder$Colour <- "black"
Tableorder$Colour[Tableorder$MeanAmpl<=1] <- "yellow"
Tableorder$Colour[Tableorder$MeanAmpl>1] <- "gold"
Tableorder$Colour[Tableorder$MeanAmpl>2] <- "orange"
Tableorder$Colour[Tableorder$MeanAmpl>3] <- "tomato"
Tableorder$Colour[Tableorder$MeanAmpl>4] <- "red"
Tableorder$Colour[Tableorder$MeanAmpl>5] <- "red2"
Tableorder$Colour[Tableorder$MeanAmpl>6] <- "red4"
Tableorder$Colour[Tableorder$MeanAmpl>7] <- "purple"
Tableorder$Colour[Tableorder$MeanAmpl>8] <- "darkorchid3"
Tableorder$Colour[Tableorder$MeanAmpl>9] <- "darkorchid4"
Tableorder$Colour[Tableorder$ActCell==0] <- "white"
col <- Tableorder$Colour

####Polygon area around all cells on the map 
chul <- as.data.frame(lapply(Tableorder[ , 2:3],"[",chull(Tableorder[ , 2:3])))
Pol <- Polygon(chul)
Areatot <- Pol@area

###Find higher Ampl cells, clustering
TableMean <- as.matrix(Table$MeanAmpl)

###'Elbow Method' to determine k for Highly active cells
fviz_nbclust(TableMean, kmeans, method = "wss")

k.max <- 10
wssAmpl <- sapply(1:k.max, 
                      function(k){kmeans(TableMean, k, nstart=50, iter.max = 15 )$tot.withinss})

plot(1:k.max, wssAmpl,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

###Apply k
TableMat <- as.matrix(Table[ , 2:3])
set.seed(Table[which(Table$MeanAmpl == max(Table$MeanAmpl)), 4])
Kmeans <- kmeans(Table[ ,1], 3)
Kmeanscluster <- as.matrix(Kmeans$cluster)
TableMatclust <- cbind(TableMat, Kmeans$cluster)
Tabledfclust <- data.frame(TableMatclust)
CellMax <- Tabledfclust[which(Table$MeanAmpl == max(Table$MeanAmpl)), 3]
TabledfclustMax <- data.frame(Tabledfclust[which(Tabledfclust$V3 == CellMax), ])

###Distance between Highly active cells
combclustMax <- combn(nrow(TabledfclustMax), 2)
connectionsclustMax <- data.frame(
  from = TabledfclustMax[combclustMax[1, ], 1:2],
  to = TabledfclustMax[combclustMax[2, ], 1:2])
connectionsclustMax$dist <- pointDistance(connectionsclustMax[ , 1:2], connectionsclustMax[ , 3:4], lonlat=F)

###Amplitude of Highly active cells
Table_MeanAmpl <- Table$MeanAmpl
ptsMax <- as.numeric(gsub("Mean","", rownames(TabledfclustMax)))
MaxCells <- (length(ptsMax) / ncol(DFfilt)) * 100
AmplMax <- as.matrix(Table_MeanAmpl[ptsMax])

###Voronoi
vtess <- deldir(x = TabledfclustMax$X, y = TabledfclustMax$Y)
Vtessdf <- as.data.frame(vtess$summary)
Vtessdf <- Vtessdf[order(Vtessdf$dir.area), ]

###Cluster, keep k = 1
set.seed(Vtessdf[which(Vtessdf$dir.area == min(Vtessdf$dir.area)),8])
kmeansarea <- kmeans(Vtessdf[ ,8], 1)
Vtessdfcluster <- cbind(Vtessdf, kmeansarea$cluster)
Clusterarea <- Vtessdfcluster[which(Vtessdfcluster[ ,8] == min(Vtessdfcluster[ ,8])), 10]
Vtessdfclust <- data.frame(Vtessdfcluster[which(Vtessdfcluster[ ,10] == Clusterarea), ])

combclust <- combn(nrow(Vtessdfclust), 2)
connectionsclust <- data.frame(
  from = Vtessdfclust[combclust[1, ], 1:2],
  to = Vtessdfclust[combclust[2, ], 1:2])
connectionsclust$dist <- pointDistance(connectionsclust[ ,1:2], connectionsclust[ ,3:4], lonlat=F)

###'Elbow Method' to determine k for nearest Highly active cells
dist <- as.matrix(connectionsclust$dist)
fviz_nbclust(dist, kmeans, method = "wss")
wssdist <- sapply(1:k.max, 
                     function(k){kmeans(dist, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wssdist,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
###Apply k
set.seed(connectionsclust[which(connectionsclust$dist == min(connectionsclust$dist)), 5])
kmeansdist <- kmeans(connectionsclust[ ,5], 3)
connectionsclustdist <- cbind(connectionsclust, kmeansdist$cluster)

###Determine Highly active cells far from the cluster
Cluster_NB <- connectionsclustdist[which(connectionsclustdist[ ,5] == min(connectionsclustdist[ ,5])), 6]
connectionsgroup <- data.frame(connectionsclustdist[which(connectionsclustdist$`kmeansdist$cluster` == Cluster_NB), ])

group1 <- connectionsgroup[ , 1:2]
row.names(group1) <- NULL
colnames(group1) <- c("x", "y")
group2 <- connectionsgroup[ , 3:4]
row.names(group2) <- NULL
colnames(group2) <- c("x", "y")
clusteredges <- unique(rbind(group1, group2))
chulMax <- as.data.frame(lapply(clusteredges[ , 1:2],"[",chull(clusteredges[ , 1:2])))

par(bg = 'white')
par(mfrow=c(1,1))
plot(Tableorder$X, Tableorder$Y, xlim=c(0, 392.66), ylim = c(0, 392.66), type="p", col="black", pch=1, cex=2.5, lwd=1, ann=FALSE, axes=FALSE)
points(Tableorder$X, Tableorder$Y, col= col, pch=16, cex=2)
polygon(chul,lty=2,border="black", lwd=1.5)
points(TabledfclustMax$X, TabledfclustMax$Y, col= "green", cex=2, lwd=2)
textxy(Tableorder$X, Tableorder$Y, Tableorder$Cells, col="black", cex=0.8)
legend("bottomright", c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), fill=c("yellow", "gold", "orange", "tomato", "red", "red2", "red4", "purple", "darkorchid3", "darkorchid4"), cex=.6, ncol=2, border=NA)
par(new=T)
plot(vtess, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1)
polygon(chulMax,lty=2,border="red", lwd=3, col=rgb(1, 0, 0, 0.25))
points(Vtessdfclust[ ,1:2], col="black", cex=2.5, lwd=2.5)

####Cluster Area
Pol <- Polygon(chul)
Polchul <- Polygon(chulMax)
Area <- (Polchul@area/Pol@area) * 100

AreaClust <- Polchul@area
Areatot <- Pol@area
RAreaClust <- AreaClust/Areatot * 100

###Determine highly active cells inside the cluster
chulMaxm <- as.matrix(chulMax)
pts <- as.matrix(Table[ , 2:3])
ptsdf <- as.data.frame(pts)
ptsInside1 <- ptsdf[in.out(chulMaxm,pts), ]
ptsInside2 <- TabledfclustMax[match(chulMax$x, TabledfclustMax$X), 1:2]
ptsInside <- unique(rbind(ptsInside1, ptsInside2))

AmplInside <- as.matrix(Table[rownames(ptsInside), 1])
ptsIn <- as.numeric(gsub("Mean","", rownames(ptsInside)))
AmplInside <- as.matrix(Table_MeanAmpl[ptsIn])
AmplInside[is.na(AmplInside)] <- 0

######Determine cells outside the cluster
ptsOutside <- data.frame(ptsdf[!(row.names(ptsdf) %in% row.names(ptsInside)), ])

AmplOutside <- as.matrix(Table[rownames(ptsOutside), 1])
AmplOutside[is.na(AmplOutside)] <- 0

ptsOut <- as.numeric(gsub("Mean","", rownames(ptsOutside)))
AmplOutside <- as.matrix(Table_MeanAmpl[ptsOut])

####Peaks
XmaxAct <- as.data.frame(mapply(FUN = function(x) {x[which(x == 0)] <- NA; x}, Xmax))
Xmaxt <- as.data.frame(t(XmaxAct))
Xmaxt$Cells <- c(1:nrow(Xmaxt))
Xmaxt$X <- Resultsxy$X
Xmaxt$Y <- Resultsxy$Y
Xmaxt$ActCell <- Table$ActCell

####Order cells with Peaks
Xmaxtorder <- Xmaxt[order(match(Xmaxt$Cells,Tableorder$Cells)), ]
Xmaxtorder$color <- Tableorder$Colour
Xmaxtorder <- Xmaxtorder[Xmaxtorder$ActCell == 1, ]

##Co-active cells
#Stim1
outmax1 <- split(Xmaxtorder, f = Xmaxtorder[ ,1])
outmax1_XY = list()
for (i in 1:length(outmax1)){
  names(outmax1)[i] <- i
  outmax1_XY[[i]] <- list(assign(paste0("outmax1_XY", i), data.frame(outmax1[[i]][(nStim+2):(nStim+3)])))
  names(outmax1_XY)[i] <- i
}
outmax1m = list()
for (i in 1:length(outmax1)){
  outmax1m[[i]] <- list(assign(paste0("outmax1m", i), mean(outmax1[[i]][ ,1])))
  names(outmax1m)[i] <- i
}
outmax1m <- as.data.frame(outmax1m)
outmax1mt <- t(outmax1m)
CoActCellmax1 = list()
for (i in 1:length(outmax1_XY)){
  CoActCellmax1[[i]] <- list(assign(paste0("CoActCellmax1", i), nrow(outmax1_XY[[i]][[1]])/nrow(Table[Table$ActCell == 1, ])*100))
}
CoActCellmax1 <- as.data.frame(CoActCellmax1)
CoActCellmax1t <- t(CoActCellmax1)

#Stim2
outmax2 <- split(Xmaxtorder, f = Xmaxtorder[ ,2])
outmax2_XY = list()
for (i in 1:length(outmax2)){
  names(outmax2)[i] <- i
  outmax2_XY[[i]] <- list(assign(paste0("outmax2_XY", i), data.frame(outmax2[[i]][(nStim+2):(nStim+3)])))
  names(outmax2_XY)[i] <- i
}
outmax2m = list()
for (i in 1:length(outmax2)){
  outmax2m[[i]] <- list(assign(paste0("outmax2m", i), mean(outmax2[[i]][ ,2])))
  names(outmax2m)[i] <- i
}
outmax2m <- as.data.frame(outmax2m)
outmax2mt <- t(outmax2m)
CoActCellmax2 = list()
for (i in 1:length(outmax2_XY)){
  CoActCellmax2[[i]] <- list(assign(paste0("CoActCellmax2", i), nrow(outmax2_XY[[i]][[1]])/nrow(Table[Table$ActCell == 1, ])*100))
}
CoActCellmax2 <- as.data.frame(CoActCellmax2)
CoActCellmax2t <- t(CoActCellmax2)

#Stim3
outmax3 <- split(Xmaxtorder, f = Xmaxtorder[ ,3])
outmax3_XY = list()
for (i in 1:length(outmax3)){
  names(outmax3)[i] <- i
  outmax3_XY[[i]] <- list(assign(paste0("outmax3_XY", i), data.frame(outmax3[[i]][(nStim+2):(nStim+3)])))
  names(outmax3_XY)[i] <- i
}
outmax3m = list()
for (i in 1:length(outmax3)){
  outmax3m[[i]] <- list(assign(paste0("outmax3m", i), mean(outmax3[[i]][ ,3])))
  names(outmax3m)[i] <- i
}
outmax3m <- as.data.frame(outmax3m)
outmax3mt <- t(outmax3m)
CoActCellmax3 = list()
for (i in 1:length(outmax3_XY)){
  CoActCellmax3[[i]] <- list(assign(paste0("CoActCellmax3", i), nrow(outmax3_XY[[i]][[1]])/nrow(Table[Table$ActCell == 1, ])*100))
}
CoActCellmax3 <- as.data.frame(CoActCellmax3)
CoActCellmax3t <- t(CoActCellmax3)

#Stim4
outmax4 <- split(Xmaxtorder, f = Xmaxtorder[ ,4])
outmax4_XY = list()
for (i in 1:length(outmax4)){
  names(outmax4)[i] <- i
  outmax4_XY[[i]] <- list(assign(paste0("outmax4_XY", i), data.frame(outmax4[[i]][(nStim+2):(nStim+3)])))
  names(outmax4_XY)[i] <- i
}
outmax4m = list()
for (i in 1:length(outmax4)){
  outmax4m[[i]] <- list(assign(paste0("outmax4m", i), mean(outmax4[[i]][ ,4])))
  names(outmax4m)[i] <- i
}
outmax4m <- as.data.frame(outmax4m)
outmax4mt <- t(outmax4m)
CoActCellmax4 = list()
for (i in 1:length(outmax4_XY)){
  CoActCellmax4[[i]] <- list(assign(paste0("CoActCellmax4", i), nrow(outmax4_XY[[i]][[1]])/nrow(Table[Table$ActCell == 1, ])*100))
}
CoActCellmax4 <- as.data.frame(CoActCellmax4)
CoActCellmax4t <- t(CoActCellmax4)

#Stim5
outmax5 <- split(Xmaxtorder, f = Xmaxtorder[ ,5])
outmax5_XY = list()
for (i in 1:length(outmax5)){
  names(outmax5)[i] <- i
  outmax5_XY[[i]] <- list(assign(paste0("outmax5_XY", i), data.frame(outmax5[[i]][(nStim+2):(nStim+3)])))
  names(outmax5_XY)[i] <- i
}
outmax5m = list()
for (i in 1:length(outmax5)){
  outmax5m[[i]] <- list(assign(paste0("outmax5m", i), mean(outmax5[[i]][ ,5])))
  names(outmax5m)[i] <- i
}
outmax5m <- as.data.frame(outmax5m)
outmax5mt <- t(outmax5m)
CoActCellmax5 = list()
for (i in 1:length(outmax5_XY)){
  CoActCellmax5[[i]] <- list(assign(paste0("CoActCellmax5", i), nrow(outmax5_XY[[i]][[1]]) / nrow(Table[Table$ActCell == 1, ]) * 100))
}
CoActCellmax5 <- as.data.frame(CoActCellmax5)
CoActCellmax5t <- t(CoActCellmax5)

outmax1min <- min(outmax1mt[ ,1])
outmax1recal <- as.data.frame(mapply(FUN=function(x) x-outmax1min, outmax1mt[ ,1]))
outmax2min <- min(outmax2mt[ ,1])
outmax2recal <- as.data.frame(mapply(FUN=function(x) x-outmax2min, outmax2mt[ ,1]))
outmax3min <- min(outmax3mt[ ,1])
outmax3recal <- as.data.frame(mapply(FUN=function(x) x-outmax3min, outmax3mt[ ,1]))
outmax4min <- min(outmax4mt[ ,1])
outmax4recal <- as.data.frame(mapply(FUN=function(x) x-outmax4min, outmax4mt[ ,1]))
outmax5min <- min(outmax5mt[ ,1])
outmax5recal <- as.data.frame(mapply(FUN=function(x) x-outmax5min, outmax5mt[ ,1]))

Xmaxtorder1 <- Xmaxtorder[order(Xmaxtorder[ ,1]), ]
Xmaxtorder1 <- Xmaxtorder1[!is.na(Xmaxtorder1$V1), ]
Xmaxtorder2 <- Xmaxtorder[order(Xmaxtorder[ ,2]), ]
Xmaxtorder2 <- Xmaxtorder2[!is.na(Xmaxtorder2$V2), ]
Xmaxtorder3 <- Xmaxtorder[order(Xmaxtorder[ ,3]), ]
Xmaxtorder3 <- Xmaxtorder3[!is.na(Xmaxtorder3$V3), ]
Xmaxtorder4 <- Xmaxtorder[order(Xmaxtorder[ ,4]), ]
Xmaxtorder4 <- Xmaxtorder4[!is.na(Xmaxtorder4$V4), ]
Xmaxtorder5 <- Xmaxtorder[order(Xmaxtorder[ ,5]), ]
Xmaxtorder5 <- Xmaxtorder5[!is.na(Xmaxtorder5$V5), ]

for (i in 1:nStim) {
  assign(paste0("XmaxtorderM", i), max(Xmaxtorder[ ,i], na.rm=T))
}
for (i in 1:nStim) {
  assign(paste0("Xmaxtorderm", i), min(Xmaxtorder[ ,i], na.rm=T))
}

nbgroupmax <- mean(c(length(outmax1),length(outmax2),length(outmax3),length(outmax4),length(outmax5)))
MaxCoActStmax <- max(c(CoActCellmax1t,CoActCellmax2t,CoActCellmax3t,CoActCellmax4t,CoActCellmax5t))

Peaks1st <- as.data.frame(outmax1recal[ ,1])
colnames(Peaks1st)[1] <-'Time (msec)'
Peaks1st$'CoActive Cells (%)'<- CoActCellmax1t[ ,1]

Peaks2nd <- as.data.frame(outmax2recal[ ,1])
colnames(Peaks2nd)[1] <-'Time (msec)'
Peaks2nd$'CoActive Cells (%)'<- CoActCellmax2t[ ,1]

Peaks3rd <- as.data.frame(outmax3recal[ ,1])
colnames(Peaks3rd)[1] <-'Time (msec)'
Peaks3rd$'CoActive Cells (%)'<- CoActCellmax3t[ ,1]

Peaks4th <- as.data.frame(outmax4recal[,1])
colnames(Peaks4th)[1] <-'Time (msec)'
Peaks4th$'CoActive Cells (%)'<- CoActCellmax4t[,1]

Peaks5th <- as.data.frame(outmax5recal[ ,1])
colnames(Peaks5th)[1] <-'Time (msec)'
Peaks5th$'CoActive Cells (%)' <- CoActCellmax5t[ ,1]

Peaks5stims <- list(Peaks1st, Peaks2nd, Peaks3rd, Peaks4th, Peaks5th) 
Peaks5stimsM <- coredata(do.call(cbind, lapply(Peaks5stims, zoo))) 
colnames(Peaks5stimsM)[c(1, 3, 5, 7, 9)] <-'Time (msec)'
colnames(Peaks5stimsM)[c(2, 4, 6, 8, 10)] <-'CoActive Cells (%)'
Peaks5stimsdf <- data.frame(Peaks5stimsM)

CoActCS1 <- ave(Peaks5stimsdf[ ,2], cumsum(Peaks5stimsdf[ ,1] == 0), FUN = cumsum)
CoActCS2 <- ave(Peaks5stimsdf[ ,4], cumsum(Peaks5stimsdf[ ,3] == 0), FUN = cumsum)
CoActCS3 <- ave(Peaks5stimsdf[ ,6], cumsum(Peaks5stimsdf[ ,5] == 0), FUN = cumsum)
CoActCS4 <- ave(Peaks5stimsdf[ ,8], cumsum(Peaks5stimsdf[ ,7] == 0), FUN = cumsum)
CoActCS5 <- ave(Peaks5stimsdf[ ,10], cumsum(Peaks5stimsdf[ ,9] == 0), FUN = cumsum)

gbdf1 <- data.frame(Peaks5stimsdf[ ,1])
colnames(gbdf1)[1] <- "Time"
gbdf1$CoAct <- CoActCS1
gbdf2 <- data.frame(Peaks5stimsdf[ ,3])
colnames(gbdf2)[1] <- "Time"
gbdf2$CoAct <- CoActCS2
gbdf3 <- data.frame(Peaks5stimsdf[ ,5])
colnames(gbdf3)[1] <- "Time"
gbdf3$CoAct <- CoActCS3
gbdf4 <- data.frame(Peaks5stimsdf[ ,7])
colnames(gbdf4)[1] <- "Time"
gbdf4$CoAct <- CoActCS4
gbdf5 <- data.frame(Peaks5stimsdf[ ,9])
colnames(gbdf5)[1] <- "Time"
gbdf5$CoAct <- CoActCS5

AvFit <- list(gbdf1, gbdf2, gbdf3, gbdf4, gbdf5)
AvFit_df <- do.call(rbind, AvFit)
AvFit_df <- na.omit(AvFit_df)
Appdf <- data.frame(approx(AvFit_df$Time, AvFit_df$CoAct, ties = mean))

###Fit NonLinear regression
fit=nls(y~y0 - a1 * exp(-b1 * x), data=Appdf, start=c(y0=99, a1=100, b1=0.016))
plot(y~x, Appdf, ylim=c(0, 100))
lines(Appdf$x, predict(fit), col='red')

summary(fit)

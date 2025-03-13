#Copyright (C) 2022-2026 F. Appaix, Univ. Grenoble Alpes, Grenoble, France.

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
library(plotrix)
library(RColorBrewer)
library(calibrate)
library(deldir)
library(gplots)
library(dplyr)
library(raster)
library(memisc)
library(mgcv)


Results <- read.table("/path/to/Data/Results.txt", header = TRUE, sep = "\t")
Resultsxy <- read.table("/path/to/Data/Resultsxy.txt", header = TRUE, sep = "\t")

###Time for 15.3Hz imaging
Fq <- 15.3
dt <- (1 / Fq)
Time <- (Results$X - 1) * dt * 1000

####Remove X column
Results$X <- NULL
Resultsxy$X.1 <- NULL

###FOR Results data Frame with NaN only!!!
Results <- Results[ ,!sapply(Results, function(x) all(is.nan(x)))]
Resultsxy <- na.omit(Resultsxy)

####Window size for rolling window and n parameter for Savitzky Golay filter 
windowSizeF0 <- Fq * 10#10 sec

#Plot a cell to define the Nb of stimulations
plot(Results$Mean30, type='l')#cell '30' for instance
#Total Nb of Stims
Nb = 7

WindowStim <- length(Time) / Nb

FindPeaksRes <- as.matrix(rollapply(Results, width = WindowStim, FUN = function(x) max(x), by = WindowStim))

row <- FindPeaksRes
for (i in 1:ncol(Results))
{
  for (j in 1:nrow(FindPeaksRes)){
    z <- which(Results[ ,i] == FindPeaksRes[j, i])
    if (length(z) > 1){
      z <- z[[1]]
      row[j, i] <- z
    }else{row[j, i] <- which(Results[ ,i] == FindPeaksRes[j, i])}
  }
}

n <- data.frame(matrix(row[1, ] - (WindowStim / 2)))
for (i in 2:(Nb)){
  n[ ,i] <- (row[i, ] - (WindowStim / 2))
}
colnames(n)[1] <- 'V1' 
for (i in 1:ncol(n))
{
  for (j in 1:nrow(n)){
    n[j, i][which(n[j, i] < 0)] <- 0
  }}

#ATTENTION n is the delay from the start of the recording to have the window get the peak properly, it can be 0!
for (i in 1:(Nb)){
  abline(v = n[30, i], col = 'red')#Row nb is the cell nb from 'plot' (line60)
}

####Create a 10sec rolling window with 50% of lowest values
F0 <- as.data.frame(rollapply(Results, width = windowSizeF0, FUN = function(z) mean(z[z < quantile(z, 0.5)]), partial = TRUE, fill = NULL, by.column = TRUE))

###DF/F0
DF <- (Results - F0) / F0
DF <- as.data.frame(mapply(FUN = function(x) ifelse(!is.na(x), x ,0), DF))

###Filtered signal
windowSizeSG <- 7
DFfilt <- as.data.frame(mapply(FUN = function(v) sgolayfilt(v, p = 3, n = windowSizeSG), DF))

####Find 5 peaks(MaxAmpl) : Y=Amplmax & X=Peaks(Index); Threshold: 2*SD, from a part of the non-filtered signal
WindowStimDF <- data.frame(matrix(row[2, ] - row[1, ]))
for (i in 2:(Nb - 1)){
  WindowStimDF[ ,i] <- (row[(i + 1), ] - row[i, ])
}
colnames(WindowStimDF)[1] <- 'V1'

nStim <- 5 ###Nb of stims kept
###Check Peaks below, if there are not 5 rows, try with a value!!
WindowStimDF1 <- mean(WindowStimDF[ ,1], na.rm = T)
#WindowStimDF1 <- 130
row7 <- mean(row[[7]], na.rm = T)

if(Nb >= 7){
  FindPeaks <- as.matrix(rollapply(DFfilt[(WindowStimDF1:row7), ], width=WindowStim, FUN = function(x) max(x), by = WindowStim))
  thres <- as.matrix(mapply(FUN = function (x) 2 * sd(x), DF[(WindowStimDF1:row7), ]))
}else{
  FindPeaks <- as.matrix(rollapply(DFfilt[(WindowStimDF1:length(Time)), ], width=WindowStim, FUN = function(x) max(x), by = WindowStim))
  thres <- as.matrix(mapply(FUN = function (x) 2 * sd(x), DF[(WindowStimDF1:length(Time)), ]))
}

Peaks <- FindPeaks
for (i in 1:ncol(DFfilt))
{
  for (j in 1:nrow(FindPeaks)){
    z <- FindPeaks[j, i]
    if(z >= thres[[i]]){
      Peaks[j, i] <- as.matrix(which(DFfilt[ ,i] == z))
    }else{
      FindPeaks[j, i] <- 0
      Peaks[j, i] <- NA
    }
  }
}

for (i in 1:ncol(DFfilt))
{
  for (j in 1:(nrow(Peaks) - 1)){
    Peaks[j, i][which(Peaks[j, i] < WindowStimDF1)] <- NA
    Peaks[(j + 1), i][which(Peaks[(j + 1), i] < (Peaks[j, i] + .5 * WindowStim))] <- NA
  }
}

for (i in 1:ncol(DFfilt))
{
  for (j in 1:nrow(FindPeaks)){
    FindPeaks[j, i][which(is.na(Peaks[j, i]))] <- NA
  }
}
Ymax <- as.data.frame(FindPeaks)
Ymax <- as.data.frame(mapply(FUN = function(x) ifelse(!is.na(x), x, 0), Ymax))

Peaksdf <- as.data.frame(Peaks)

####Find Time (Xmax) for each peak
Xmax <- data.frame()
for (i in 1:ncol(DFfilt))
{
  for (j in 1:nrow(Peaks)){
    if (!is.na(Peaksdf[j, i])){
      Xmax[j, i] <- lapply(FUN = function(x) Time[x], Peaksdf[j, i])
    }else{
      Xmax[j, i] <- 0}
  }
}
colnames(Xmax) <- colnames(DFfilt)

####Find peaks start Y=Yst & X=Starts(Index)
DFfiltmean <- data.frame()
for (i in 1:ncol(DFfilt)){
  for (j in 1:Nb){
    DFfiltmean[j, i] <- as.data.frame(mean(DFfilt[(n[i, j]:(n[i, j] + (windowSizeF0 / 3))), i]))
  }
}
colnames(DFfiltmean) <- colnames(DFfilt)

Starts <- Peaksdf
for (i in 1:ncol(DFfilt))
{
  for (k in 1:nrow(Peaksdf)){
    if((Peaksdf[k, i] >= mean(Peaks[1, ], na.rm = T) / 2) & (!is.na(Peaksdf[k, i]))){
      z <- Peaksdf[k, i]
      Startstemp <- DFfilt[z, i]
      while (Startstemp > DFfiltmean[k + 1, i])
      {
        z <- z - 1
        Startstemp <- DFfilt[z, i]
      }
      val20hz <- DFfilt[(z - 1), i]
      while (val20hz <= Startstemp)
      {
        z <- z - 1
        Startstemp <- DFfilt[z, i]
        val20hz <- DFfilt[(z - 1), i]
      }
      Starts[k, i] <- z
    }else{
      Starts[k, i] <- 1}
  }}

Startsdf <- as.data.frame(Starts)

Yst <- Startsdf
for (i in 1:ncol(Startsdf)){
  for (j in 1:nrow(Startsdf)){
    if (Startsdf[j, i] != 1){
      Yst[j, i] = DFfilt[(Startsdf[j, i]), i]
    }else{
      Yst[j, i] <- 0}
  }}
names(Yst) <- names(Startsdf)

####Find Time (Xst) for each peak start  
Xst <- as.data.frame(lapply(FUN = function(x) Time[x], Startsdf))

## Mean ampl MSNs DMS all animals 
MeanAmpl_MSN <- 3.3535#Calculated value

###Mean Ampl Peak/cell
Ydelta <- as.data.frame(mapply('-', Ymax, Yst))
Ydelta[Ydelta == 0] <- NA
Mean <- as.matrix(apply(Ydelta, 2, mean, na.rm = TRUE))
SD <- as.matrix(sapply(FUN = function(x) std.error(x), Ydelta))

###Mean RiseTime Peak/cell
Xdelta <- as.data.frame(mapply('-', Xmax, Xst))
Xdelta[Xdelta == 0] <- NA
Meanrise <- as.matrix(apply(Xdelta, 2, mean, na.rm = TRUE))
SDrise <- as.matrix(sapply(FUN = function(x) std.error(x), Xdelta))

####Find End of Peaks
Ends <- Peaksdf
for (i in 1:ncol(DFfilt))
{
  for (k in 1:nrow(Peaksdf)){
    z <- Peaksdf[k, i]
    if(!is.na(Peaksdf[k, i])){
      Endstemp <- DFfilt[z, i]
      while (Endstemp >= DFfiltmean[(k + 1), i])
      {
        z <- z + 1
        Endstemp <- DFfilt[z, i]
      }
      val20hz <- DFfilt[(z + 1), i]
      while (val20hz <= Endstemp)
      {
        z <- z + 1
        Endstemp <- DFfilt[z, i]
        val20hz <- DFfilt[(z + 1), i]
      }
      Ends[k, i] <- z
    }}}

Endsdf <- as.data.frame(Ends)

Yend <- Endsdf
for (i in 1:ncol(Endsdf)){
  for (j in 1:nrow(Endsdf)){
    if (Endsdf[j, i] != 1){
      Yend[j, i] = DFfilt[(Endsdf[j, i]), i]
    }else{Yend[j, i] <- 0}
  }}
names(Yend) = names(Endsdf)

####Find Time (XEnd) for each peak end
Xend <- as.data.frame(lapply(FUN = function(x) Time[x], Endsdf))

Min <- mapply(FUN = function (x) min(x), DF)
Min <- min(Min)
max <- mapply(FUN = function (x) max(x), DF)
max <- max(max)
if(Nb >= 7){
  th <- as.matrix(mapply(FUN = function (x) (mean(x) + 2 * sd(x)), DF[c(WindowStim:((nStim + 1) * WindowStim)), ]))#keep the same window as above, Findpeaks part.
}else{
  th <- as.matrix(mapply(FUN = function (x) (mean(x) + 2 * sd(x)), DF[c(WindowStim:(Nb * WindowStim)), ]))#keep the same window as above, Findpeaks part.
}


###Interneurons = indicate the #cells of interneurons
In <- matrix(c(1, 2))#N° of the columns, can be several cells

###Time For the Decay, Time for the peak response
DTime <- as.data.frame(mapply('-', Xend, Xmax))
DTime[DTime == 0] <- NA
DecayTime <- as.matrix(apply(DTime, 2, mean, na.rm = T))
SDDecayTime <- as.matrix(sapply(FUN = function(x) std.error(x), DTime))

####Create 'Table' with (X, Y), AmplPeak (Mean, SD), RiseTime (Mean, SD) and Cells
MeanAmpl <- as.data.frame(Mean, row.names = NULL)
Table <- MeanAmpl
colnames(Table)[1] <-'MeanAmpl'
Table$SD <- SD
Table$RiseTime <- Meanrise
Table$SDrise <- SDrise
Table$DecayTime <- DecayTime
Table$SDDecayTime <- SDDecayTime
Table$X <- Resultsxy$X
Table$Y <- Resultsxy$Y
Table$Cells <- c(1:nrow(Mean))

####Add column 'Cell types' in Table
Table$Celltype <- paste0("MSN", Table[c(1:nrow(Mean)), 10])
Table[In, 10] <- "InterNeu"
Table[which(Table$Celltype == "InterNeu"), 1] <- 0

###Identify the SST IN on the map and add SST calculated values
SST <- matrix(c(1, 2))#N° for SST
SSTdf <- data.frame(Table[SST, 7:9])
Table[SST, 10] <- "SST"

####Add column 'Active Cells' in Table
th <- round(th, 1)
Meanr <- round(Mean, 1)
pm <- as.matrix(Meanr > th)
Table$ActCell <- as.numeric(pm)

####Map Active cells vs Total detected cells
ActCell <- which(Table$ActCell == 1)
NonActCell <- which(Table$ActCell == 0)
Table[NonActCell, 1] <- 0
Table[which(Table$Celltype == "InterNeu"), 11] <- 0
RActCell <- (length(ActCell) / length(Resultsxy$X)) * 100 # % act cells

###Order Table with Ampl
Tableorder <- Table[order(Table[ ,1]), ]

###Create colorAmpl palette
Tableorder$Colour <- "black"
Tableorder$Colour[Tableorder$MeanAmpl <= 1] <- "yellow"
Tableorder$Colour[Tableorder$MeanAmpl > 1] <- "gold"
Tableorder$Colour[Tableorder$MeanAmpl > 2] <- "orange"
Tableorder$Colour[Tableorder$MeanAmpl > 3] <- "tomato"
Tableorder$Colour[Tableorder$MeanAmpl > 4] <- "red"
Tableorder$Colour[Tableorder$MeanAmpl > 5] <- "red2"
Tableorder$Colour[Tableorder$MeanAmpl > 6] <- "red4"
Tableorder$Colour[Tableorder$MeanAmpl > 7] <- "purple"
Tableorder$Colour[Tableorder$MeanAmpl > 8] <- "darkorchid3"
Tableorder$Colour[Tableorder$MeanAmpl > 9] <- "darkorchid4"
Tableorder$Colour[Tableorder$ActCell == 0] <- "white"
col20hz <- Tableorder$Colour

####Polygon area from cells map 
chul20hz <- as.data.frame(lapply(Tableorder[ ,7:8],"[",chull(Tableorder[ ,7:8])))
Pol20hz <- Polygon(chul20hz)
Areatot20hz <- Pol20hz@area

###Find HA cells using threshold method
TableHA <- Tableorder
TableHA[ ,1][TableHA[ ,1] < MeanAmpl_MSN] = 0
TableHA <- TableHA[which(TableHA[ ,1] != 0), ]
TableHAlxy <- data.frame(TableHA[ ,7:8])

###Amplitude of HA cells
Table_MeanAmpl <- Table$MeanAmpl
pts_HA <- as.numeric(gsub("Mean", "", rownames(TableHA)))
Ampl_HA <- as.matrix(Table_MeanAmpl[pts_HA])
Mean_Ampl_HA <- mean(Ampl_HA)

###LA cells
pts_LA <- Table[!(Table$Cells %in% pts_HA), 9]
Ampl_LA <- as.matrix(Table_MeanAmpl[pts_LA])
Mean_Ampl_LA <- mean(Ampl_LA)


###MeanAmpl MSN within the slice
MeanAmpl_MSN <- mean(Tableorder[which(Tableorder[ ,11] != 0), 1])

###SST DF -> DFfilt
Fs <- 15.3
bf1 <- butter(4, 1/(Fs/2), type="low")
SST_DF <- DF[ , SST]
SST_DFfilt <- filtfilt(bf1, SST_DF)

###SST find start and Peak
par(mfrow=c(1,1))
SST_DFfilt <- as.matrix(SST_DFfilt)
Nb = 7
WindowStim <- 1000/Nb

FindPeaks_SST_DFfilt <- rollapply(SST_DFfilt, width=WindowStim, FUN=function(x) max(x), by=WindowStim)
Peaks_SST_DFfilt <- FindPeaks_SST_DFfilt
for (j in 1:nrow(FindPeaks_SST_DFfilt)){
  z <- which(SST_DFfilt[] == FindPeaks_SST_DFfilt[j])
  if (length(z) > 1){
    z <- z[[1]]
    Peaks_SST_DFfilt[j] <- z
  }else{Peaks_SST_DFfilt[j] <- which(SST_DFfilt[ ] == FindPeaks_SST_DFfilt[j])}
}

Wind_SST_DFfilt <- data.frame(matrix(Peaks_SST_DFfilt[1] - (WindowStim/2)))
for (i in 2:(Nb)){
  Wind_SST_DFfilt[i] <- (Peaks_SST_DFfilt[i] - (WindowStim/2))
}
colnames(Wind_SST_DFfilt)[1] <- 'V1' 

Findstart_SST_DFfilt <- FindPeaks_SST_DFfilt
###SOIT
for (i in 1:5){
  Findstart_SST_DFfilt[i, ] <- min(SST_DFfilt[(Wind_SST_DFfilt[ ,i+1]+10):(Wind_SST_DFfilt[ ,i+1]+10), ])
}
###OU le plus simple, les callées entre les pics et on prend le 1er
Findstart_SST_DFfilt[1,] <- min(SST_DFfilt[0:Peaks_SST_DFfilt[1,],])
for (i in 2:5){
  Findstart_SST_DFfilt[i,] <- min(SST_DFfilt[Peaks_SST_DFfilt[i-1,]:Peaks_SST_DFfilt[i,],])
}
###OU
Findstart_SST_DFfilt[1,] <- min(SST_DFfilt[(Wind_SST_DFfilt[,1]):(Wind_SST_DFfilt[,2])])
for (i in 1:5){
  Findstart_SST_DFfilt[i,] <- min(SST_DFfilt[Peaks_SST_DFfilt[i,1]:Peaks_SST_DFfilt[i+1,]])
}
###OU celui-là c'est si on ne prend pas le 1er pic
Findstart_SST_DFfilt[1,] <- min(SST_DFfilt[(Wind_SST_DFfilt[,2]+10):(Wind_SST_DFfilt[,2]+(Wind_SST_DFfilt[,3]-Wind_SST_DFfilt[,2])/2),])
for (i in 2:5){
  Findstart_SST_DFfilt[i-1,] <- min(SST_DFfilt[(Wind_SST_DFfilt[,i]+10):(Wind_SST_DFfilt[,i]+(Wind_SST_DFfilt[,i+1]-Wind_SST_DFfilt[,i])/2),])
}
Findstart_SST_DFfilt[5,] <- min(SST_DFfilt[(Wind_SST_DFfilt[,6]):Peaks_SST_DFfilt[6,],])


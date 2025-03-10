1. System requirements

Install last release of R (version 3.6.3 or 4.0.2) and Rstudio (1.3.959) that can be found here:
-R: https://www.r-project.org
-Rstudio: https://rstudio.com/products/rstudio/download/
Warning - version of R (3.6.1) from anaconda 1.9.12 is not compatible

The presented code was used on a Windows operating system (Windows 10 and 7)

Libraries corresponding to the same version of Rstudio were used

2. Installation guide

Unzip received folder containing an R file and the .txt files ('Results.txt' and 'Resultsxy.txt')
Load the 'Code.R' file in Rstudio

3. Demo

Instructions to run on data

The files needed are: 'Results.txt', 'Resultsxy.txt', 'Code.R'
Once Rstudio is open, load the 'Code.R' by 'Opening an existing file'
Install the packages necessary to run the script
  > install.packages("library_name")
Load the following packages by running the script:
  > zoo (v 1.8.8)
  > signal (v 0.7.6)
  > grDevices (v 4.0.2)
  > sp (v 1.4.2)
  > factoextra (v 1.0.7)
  > cluster (v 2.1.0)
  > fpc (v 2.2.5)
  > RColorBrewer (v 1.1.2)
  > calibrate (v 1.7.7)
  > deldir (v 0.1.25)
  > gplots (v 3.0.3)
  > dplyr (v 1.0.0)
  > raster (v 3.1.5)
  > memisc (v 0.99.22)
  > mgcv (v 1.8.31)
Adjust the link towards the downloaded .txt files to your local directory
  > example: "C:/Users/MyName/Documents/example.txt"

Expected output

The R function 'View()' can be used to visualize the different generated tables
The functional map resulting from the analysis is shown in the bottom right window
Plots of data in the code can be viewed using the 'Plot()' function

A typical install time on a 'normal' computer (between 8 and 16 Go RAM) is between 5 and 10 minutes

4.Instructions for use

The sent .txt files correspond to a set of data
  > "Results.txt" file contains the intensity of the cells in function of samples (one row one sample, one column one sample)
  > "Results.txt" file contains the coordinates of the cells
The instructions described in '3.' apply for this data set

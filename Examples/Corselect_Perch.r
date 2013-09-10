################################################################################################
# 
# Corselect.r - Simultaneous estimation of selectivity parameters for Gillnets & Cormorants 
# Application to European Perch Data from Curonian Lagoon, Lithuania 
# 
# by Athol Whitten (awhitten@gmail.com)
# Melbourne, Australia, 2012
#
################################################################################################

#Remove all objects currently in workspace as precaution;
rm(list=ls())

#Specify directory of working folder;
folder <- "C:/Corselect"

#Specify folder for printing plots of results (set here as subfolder of working folder named 'Output';
output.folder <- paste(folder,"/Output",sep="")

#Source Corselect functions;
source(paste(folder,"/Corselect.r",sep=""))

#Apply Corselect functions to Perch Data. Get the data and set the meshsizes vector;
perch.data <-(read.table(paste(folder,"/Data/Perch.dat",sep=""),header=TRUE))
perch.meshsizes <- c(5,14,17,21.5,25,30,33,38,45)

#View and check the data and meshsize vector;
print(perch.data)
print(perch.meshsizes)

#Set initial parameter values for Theta 1 through Theta 4;
gn.parms <- c(1,10) #Theta 1 and 2.
gc.parms <- c(2,0.5) #Theta 3 and 4.

#Implement optimisation with all data, using Nelder-Mead method as part of R optim function:
perch.fit <- optim(c(gn.parms,gc.parms),nllhood,cdata=perch.data,meshsizes=perch.meshsizes,gn.sel=gn.sel,gc.sel=gc.sel,hessian=T,method="Nelder-Mead")

#Get parameter estimates and the mode and variance estimates for the first two nets and the Cormorants;
estimates(fit=perch.fit,meshsizes=perch.meshsizes)

#Get plots of base data (barplot) and of the estimated selectivity curves over a specified range of plot lengths;
plots(fit=perch.fit,meshsizes=perch.meshsizes,cdata=perch.data,BS=TRUE,plot.lens=seq(0.0001,50,0.01),label=TRUE,save=TRUE,save.to=output.folder,name="perch")

#End of script.
################################################################################################

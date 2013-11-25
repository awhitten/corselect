#####################################################################################################################
#
# 	Corselect.r: R Functions for the simultaneous estimation of selectivity parameters for Gillnets & Cormorants.
#	
# 	In using these methods an assumed selectivity function (a gamma function) is fitted directly to catch 
#	data for a number of gill-nets with specific mesh-sizes, with the parameters of the selectivity function 
#	being estimated simultaneously across mesh-sizes and length-classes. An additional selectivity function 
#	(a log-normal function) is simultaneously fit to 'catch' data for Great Cormorants.
#
# 	Code by Athol Whitten (awhitten@gmail.com), 
# 	Melbourne, Australia, 2012
#
# 	This code is free to use (see LICENSE.txt), but please give credit to the authors of the study for which 
#   the code was written: See Troynikov et al, 2013: Cormorant catch concerns for fishers: Estimating the 
# 	size-selectivity of a piscivorous bird, PLOS ONE. DOI: 10.1371/journal.pone.0077518
#
#######################################################################################################################

#Define Selectivity Function (S1): Relative selectivities of gill-nets using a functional form of the standard Gamma function (with modal value rescaled to one);
gn.sel <- function(length.class,meshsize,theta){
	B=-0.5*(theta[1]*meshsize-((theta[1]^2)*(meshsize^2)+4*theta[2])^0.5)
	A=(theta[1]*meshsize)/B
	sel=((length.class/(A*B))^A)*exp(A-(length.class/B))
	return(sel)
	}

#Define Selectivity Function (S2): Relative selectivities of Cormorants using a functional form of the standard Log-Normal function (with modal value rescaled to one);	
gc.sel <- function(length.class,theta){
	sel=exp((-(log(length.class)-theta[3])^2)/(2*theta[4]^2) - (theta[4]^2)/2 + theta[3])/length.class
	return(sel)
	}

#Define a log-likelihood function to be maximised (implemented as a negative log-likelihood and minimised);
nllhood <- function(theta,cdata,meshsizes,gn.sel,gc.sel){
	catch.data=cdata[which(apply(cdata[,-1],1,sum)>0),-1]
	length.class=cdata[which(apply(cdata[,-1],1,sum)>0),1]
	gn.sels=outer(length.class,meshsizes,gn.sel,theta)
	gc.sels=gc.sel(length.class,theta)
	smatrix=cbind(gn.sels,gc.sels)
	mu=apply(catch.data,1,sum)/apply(smatrix,1,sum)
	nllhood=-sum(catch.data*(log(mu*smatrix))-mu*smatrix)
	return(nllhood)
	}	

#Create function to get estimates from particular fit:
estimates <- function(fit,meshsizes){
	theta=fit$par
	B1=-0.5*(theta[1]*meshsizes[1]-((theta[1]^2)*(meshsizes[1]^2)+4*theta[2])^0.5)
	A1=(theta[1]*meshsizes[1])/B1
	B2=-0.5*(theta[1]*meshsizes[2]-((theta[1]^2)*(meshsizes[2]^2)+4*theta[2])^0.5)
	A2=(theta[1]*meshsizes[2])/B2
	Net1.Mode=A1*B1
	Net2.Mode=A2*B2
	Net.Variance=(A1+1)*B1^2
	Net.StDev=sqrt(Net.Variance)
	GC.Mode=exp(theta[3]-(theta[4]^2))
	GC.Variance=exp(2*theta[3]+(theta[4]^2))*(exp(theta[4]^2)-1)
	parmest=rbind("Theta 1"=fit$par[1],"Theta 2"=fit$par[2],"Theta 3"=fit$par[3],"Theta 4"=fit$par[4])
	stnderr=sqrt(diag(solve(fit$hessian))) #The standard error of a parameter estimates is the square root of the diagonal of the inverse of the hessian matrix.
	output.a=cbind(parmest,stnderr)
	colnames(output.a)=c("Estimate","S.E.")
	output.b=rbind("Gillnet Mode (Mesh 1)"=Net1.Mode,"Gillnet Mode (Mesh 2)"=Net2.Mode,"Gillnet Variance"=Net.Variance,"Gillnet StDev"=Net.StDev,"Cormorant Mode"=GC.Mode,"Cormorant Variance"=GC.Variance)
	colnames(output.b)="Estimate"
	print(output.a)
	print(output.b)
	}

#Create function to get plots of data and of selectivity curves;
plots <- function(fit,meshsizes,cdata,plot.lens,label=TRUE,save=FALSE,save.to="directory",name="new",BS="TRUE"){
	theta=fit$par
	gn.sels=outer(plot.lens,meshsizes,gn.sel,theta)
	gc.sels=gc.sel(plot.lens,theta)
	smatrix=cbind(plot.lens,gn.sels,gc.sels)
	
	windows(record=TRUE,width=900,height=500)
	
	barplot(as.matrix(cdata[,-1]),beside=TRUE,ylab="Frequency",xlab="Gear Type",ylim=c(0,1.1*max(cdata[,-1])))
	
	if(BS==TRUE){
		matplot(plot.lens,gn.sels[,-1],type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
		matlines(plot.lens,gn.sels[,1],type="l",col=1,lty=3)
		}
	
	if(BS==FALSE){
		matplot(plot.lens,gn.sels,type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
		}
		
	matlines(plot.lens,gc.sels,type="l",col=1,lwd=2)
		
		if(label==TRUE){
			for(i in 1:ncol(cdata[-1])){
			text(smatrix[which(smatrix[,i+1]==max(smatrix[,i+1])),1],1.03,colnames(cdata[-1])[i],cex=0.8)
			}
		}
	
	if(save==TRUE){
			postscript(file=paste(save.to,"/",name,".data.eps",sep=""),onefile=FALSE, horizontal=TRUE, pointsize=12, bg="white")
			barplot(as.matrix(cdata[,-1]),beside=TRUE,ylab="Frequency",xlab="Gear Type",ylim=c(0,1.1*max(cdata[,-1])))
			dev.off()
			
			postscript(file=paste(save.to,"/",name,".selectivity.eps",sep=""),onefile=FALSE, horizontal=TRUE, pointsize=12,bg="white")
			
			
				if(BS==TRUE){
					matplot(plot.lens,gn.sels[,-1],type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
					matlines(plot.lens,gn.sels[,1],type="l",col=1,lty=3)
					}
	
				if(BS==FALSE){
				  matplot(plot.lens,gn.sels,type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
				  }
			
			
			matlines(plot.lens,gc.sels,type="l",col=1,lwd=2)
				
				if(label==TRUE){
					for(i in 1:ncol(cdata[-1])){
					text(smatrix[which(smatrix[,i+1]==max(smatrix[,i+1])),1],1.03,colnames(cdata[-1])[i],cex=0.8)
					}
				}
			dev.off()
			}
}

#######################################################################################################################
#End of Corselect.r.
#######################################################################################################################

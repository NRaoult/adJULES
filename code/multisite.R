#!/usr/bin/env Rscript

source(file.path(code,'fhjules.R')) # load adjules
source(file.path(code,'multifunctions.R')) # loads in multisite_functions

multisite <- 'y'
multichoice <- set
########################################################
#
# Set up multisite function
#
########################################################

nsites <- length(multichoice)
multijinfile <- rep("",nsites)
multidatfile <- rep("",nsites)
multijinname <- locs[multichoice]
print(multijinname)
multiw <- array(NA,dim=c(nts,nts,nsites))

for (i in 1:nsites){
   multijinfile[i] <- paste(verifydir,"/control/",locs[multichoice[i]],"_",exp_year[multichoice[i]],".jin",sep="") 
   multidatfile[i] <- paste(verifydir,"/data/obsts-",locs[multichoice[i]],"_",exp_year[multichoice[i]],".da",sep="")
   multiw[,,i]     <- w
}

xoffmulti   <- array(NA,dim=c(n,nsites))
for (ilocal in 1:nsites){
   unlink(  c("modts.da" , "obsts.da" , "control.jin", "avmodts.da", "avobsts.da")  )                   
   temp1 = file.copy( multidatfile[ilocal] , "obsts.da" ) ; if (!temp1) stop("Error: observations not copied over")
   temp2 = file.copy( multijinfile[ilocal],"control.jin" ); if (!temp2) stop("Error: jinfile not copied over")
   rm(temp1 , temp2)

   temp <- .Fortran("initfunc",n=n,x=double(n))
   xoff <- temp$x
   rm(temp)
   .Fortran("postfunc",n=n,m=m,x=xoff,y=double(m))      
   
   xoffmulti[,ilocal]    <- xoff
}


# starting point
zoff          <- xoff[xvary] 
names(zoff)   <- znames
nz            <- length(zoff)

if (update_start == 'on'){

	znudged <- as.numeric(z_av) # averaged starting point over PFTs
	nfactor <- 'same_start'
	check1 <- sum(znudged < zlower) ; check2 <- sum(znudged > zupper)
  if (check1+check2 > 0) print('Out of bounds')

	x1   <- xoff

	zoff <- as.numeric(znudged)
	x1[xvary] <- zoff

	for (i in 1:nsites){
		xoffmulti[parameterchoice,i] <- zoff 
		x1 <- xoffmulti[,i]
   	jinfile <- paste(verifydir,"/control/",locs[multichoice[i]],"_",exp_year[multichoice[i]],".jin",sep="")
    jinname <- locs[multichoice[i]]
		source(file.path(code, 'update_control.R'))
   	multijinfile[i] <- paste0(verifydir,"/control/optimised_",locs[multichoice[i]],".jin") 	
	}
}


# performs optimisation
source(file.path(code,'multifindopt_save.R'))



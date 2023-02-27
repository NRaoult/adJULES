# *****************************COPYRIGHT*******************************
# (c) University of Exeter 2011 . All rights reserved.
#
# contact: Tim Jupp t.e.jupp@exeter.ac.uk
#
# This routine has been licensed to the other JULES partners for use
# and distribution under the JULES collaboration agreement, subject
# to the terms and conditions set out therein.
#
# [Met Office Ref SC0237] 
# *****************************COPYRIGHT*******************************
#
# Description:
# 
# R code to perform parameter optimisation
# run this after adjules.R
# CHECK: removed x0 def in clean up

ztemp <- zoff

print(tsuse)
print(ztemp)
print(paste("Initial cost (with *initial* weight matrix) is", multif(ztemp, lgetw = TRUE))) # equivalent to getw and save w

tol      <- 1e-5 #10
savew    <- list()
savez    <- list()
savehess <- list() 
saveerr  <- list()
savezsample  <- list()
 
opt <- multiminim(ztemp)      # do a first pass on optimisation, this calls multif(z, lgetw = FALSE)

j <- 0 
jmax <-50
if (ss == 26) jmax <- 5

while (abs(sum(opt$par-ztemp))> tol && j < jmax){         # iterate a few times, finding minimum and resetting w
	        print(opt$counts)
		j <- j+1
	        print(paste('Pass number',j,'complete'))     
	        savew[[j]]    <- multiw 
	        savez[[j]]    <- opt$par
		z1 <- opt$par 
	        saveerr[[j]]  <- err[,c(3,5)]      # probably only works for single  
	        ztemp         <- opt$par            # store temporary optimum
	        multif(ztemp, lgetw = TRUE)       # ensure optimal time-series saved to disk and update w matriz
		   # savehess[[j]] <- symmpart(multid2f(ztemp))

			#hess     <- symmpart(
			 #           diag(1./(zupper-zlower)) %*% 
			  #          nearPD(   diag(zupper-zlower) %*% savehess[[j]] %*% diag(zupper-zlower)    )$mat %*%
			   #         diag(1./(zupper-zlower)) )

				#nsample <- 10000
				#zsample <- rtmvnorm(n=nsample, 
			         #            mean=z1, 
				#	                H=hess,
			         #           lower=zlower,
				#	            upper=zupper,
				#	        algorithm="gibbs",
				#	      start.value=z1
					         #      )
			#savezsample[[j]] <- zsample
	        opt <- multiminim(ztemp)

}

#########################################################################
#
# Save results
#
#########################################################################

convergence <- opt$convergence   # store diagnostic info
message     <- opt$message       # about optimisation
hess_num    <- opt$hessian       # store numerical hessian

z1          <- opt$par        # store final optimum as a z-vector...
x1          <- xoff       #
x1[xvary]   <- z1             # ...and as an x-vector

#
# now store differences between z0 and z1 runs
#

z0 <- zoff

params <- cbind(seq(1,nz),z0,z1,rep(NA,nz),z1==zlower,z1==zupper) # see old and new parameters side by side

rownames(params) <- xnames[xvary]
colnames(params) <- c("zw","initial","final","perc","z1==zlower","z1==zupper")
params[,4] <- 100 * (params[,3] - params[,2])/(zupper-zlower)

print(zapsmall(params))
print(opt$message)
print(opt$convergence)
#print(paste("With **final** weight matrix, cost is reduced from",multif0,"to",multif1))
hess_taf <- symmpart(multid2f(z1))

hess     <- symmpart(
            diag(1./(zupper-zlower)) %*% 
            nearPD(   diag(zupper-zlower) %*% hess_taf %*% diag(zupper-zlower)    )$mat %*%
           diag(1./(zupper-zlower)) )

	nsample <- 10000
	zsample <- rtmvnorm(n=nsample, 
                     mean=z1, 
		                H=hess,
                    lower=zlower,
		            upper=zupper,
		        algorithm="gibbs",
		      start.value=z1
		               )
	z1lo <- NA * z1
	z1hi <- NA * z1	      
	      
	for (i in 1:nz){
   		z1lo[i] <- quantile(zsample[,i], prob=0.1)
   		z1hi[i] <- quantile(zsample[,i], prob=0.9)
	}

parmat <- (cbind(zlower,z1lo,z0,z1,z1hi,zupper))
colnames(parmat) <- c("lower","quantile_0.1","original","optimal","quantile_0.9","upper")

#if (!(type %in% c('single','ran5','double','chosen'))) ss <- ''
#if (update_start == 'on') ns <- nfactor
#save.image(paste0(result,'zsamples_8_',experiment,'_',type,ss,'_',nfactor,'_',alpha,'_',stream,'.RData'))
#dg <- 'off' if (diagonal == 'TRUE') dg <- 'on'

if (update_start == 'on'){
	us <- '_samestart'
} else {
	us <-''
}
save.image(paste0(result,'zsamples_8_',experiment,'_',type,ss,'_canopy_',canopy,'_',alpha,'_',stream,'_wsave_',windtype,us,'.RData'))

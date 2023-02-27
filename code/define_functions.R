load_site <- function(i = 1){ #or setup_site? what the best name for this script?

  jinname <- locs[i]
  print(jinname)
  jinfile <- paste(verifydir,"/control/",locs[i],"_",exp_year[i],".jin",sep="")
  datfile <- paste(verifydir,"/data/obsts-",locs[i],"_",exp_year[i],".da",sep="") #,

  unlink(  c("modts.da" , "obsts.da" , "control.jin", "avmodts.da", "avobsts.da")  )                   # DP change
    temp1 = file.copy( datfile , "obsts.da" ) ;  if (!temp1) stop("Error: observations not copied over")
    temp2 = file.copy( jinfile,"control.jin" ) ; if (!temp2) stop("Error: jinfile not copied over")
  rm(temp1 , temp2)

}

change_canopy <- function(canopy){
control_file     <- scan('control.jin',character(0),sep="\n",comment.char="#")
new_control_file <- control_file # create new control file object

## change canopy model
can_line <- grep('can_rad',control_file)
new_control_file[can_line] <- sub(1,canopy,control_file[can_line])
write(new_control_file,"control.jin")
print(paste('control.jin file ammended into canopy model',canopy))
}

read_timeseries <- function(bfile='obsts.da', # binary file to be read in
			     nmax= as.integer(1e7), # maximum number of data points to read in - could change if required but may increase storage requirement
				 			 nts = 6) 
{
	con  <- file(bfile,"rb") #read in binary file
	temp <- readBin(con =con,
	                what=numeric(),
		            size=8,
			      endian="little",
		               n=nmax)
	if (length(temp)==nmax){print("readts WARNING may have hit maximum data size")}
	close(con)
	nt2 <- length(temp)/nts

	timeseries <- array(NA,dim=c(nt2,nts))
	for (i in 1:nts){
    	timeseries[,i] <- temp[seq(i,nt2*nts-nts+i,nts)]
	}
	rm(temp)
	timeseries[timeseries == -9999.99] <- NA   # reset NA from dble-valued NA flag
    out = timeseries
    out
}

find_modts <- function(z=z,xoff=xoff){
  f(z,xoff)
  out <- read_timeseries('avmodts.da')
  out
}


cal_rmse <- function(err=err,its=1){ 
  out=sqrt(mean((err[,its])^2,na.rm=TRUE))
  out
}


###
# functions of measure
###
find_modts <- function(z=z,xoff=xoff){
  system("rm avmodts.da")
  print(f(z,xoff))
  out <- read_timeseries('avmodts.da')
  out
}

normalised_rmse <- function(modts=modts,its=3){
  RMSED <- cal_rmse((modts-obsts),its)/nrow(modts-obsts)
  ymin <- min(sapply(modts_set,function(modts) min(modts[,its])))
  ymax <- max(sapply(modts_set,function(modts) max(modts[,its])))
  out  <- RMSED/(ymax-ymin)
  out
}

FVU <- function(modts=modts,its=3){
  o_bar <- mean(obsts[,its],na.rm=TRUE)
  top <- sum((obsts[,its]-modts[,its])^2,na.rm=TRUE)
  bottom <- sum((obsts[,its]-o_bar)^2,na.rm=TRUE)
  out  <- top/bottom
}

measure_rmse <- function(modts=mosts,dp=5){
  LE_RMSED <- normalised_rmse(modts,3)
  GPP_RMSED <- normalised_rmse(modts,5)
  out <- round((LE_RMSED+GPP_RMSED)/sqrt(2),dp)
  out
}

measure <- function(modts=mosts,dp=5){
  LE_RMSED <- FVU(modts,3)
  GPP_RMSED <- FVU(modts,5)
  out <- round((LE_RMSED+GPP_RMSED)/2,dp)
  out
}

rmse <- function(modts=modts,its=3){
err <- modts[,its]-obsts[,its]
out=sqrt(mean((err)^2, na.rm=TRUE))
out
}
rmseGPP <- function(modts=modts,its=5){
err <- modts[,its]-obsts[,its]
out=sqrt(mean((err)^2, na.rm=TRUE))
out
}

measure2 <- function(modts,dp=5){
RMSED_LE_orig <- cal_rmse((modts0-obsts),3)/nrow(modts0-obsts)
RMSED_GPP_orig <- cal_rmse((modts0-obsts),5)/nrow(modts0-obsts)
RMSED_LE_new <- cal_rmse((modts-obsts),3)/nrow(modts-obsts)
RMSED_GPP_new <- cal_rmse((modts-obsts),5)/nrow(modts-obsts)
LE <- 1- (RMSED_LE_orig-RMSED_LE_new)/RMSED_LE_orig
GPP <- 1- (RMSED_GPP_orig-RMSED_GPP_new)/RMSED_GPP_orig
out <- round((LE+GPP)/2,dp)
out 
}


find_parameters_to_remove <- function(hessian){#lambda, pft,type='',datastream=''){ #data = _NEE for NEE experiments
  #/scratch/nr278/Dropbox/finding_lambda.R
  #load(paste0('Results/zsamples_all',pft,'_',type,lambda,datastream,'.RData'))
  if (is.positive.definite(hessian)){
    print('already positive definite')
    working_sets <- NA
  } else {
    c <- 0
    parameters <- seq(1,8,1)
    working_sets <- list(NA)
    for (i in 1:7){
      subset <- combn(parameters,i)
        for (s in 1:ncol(subset)){
          if (is.positive.definite(hessian[-subset[,s],-subset[,s]])){
            c <- c+1
            working_sets[[c]] <- subset[,s]
          }
        }
        if (!is.na(working_sets)[1]) break  
    }
  }
  out <- working_sets
  out
}





#################
#
# script defining functions to generate plots after optimisations
#
#################
numdf1 <- function(z,eps,i)  # 1st order finite difference approximation to df at z
{      e <- 0 * z
       e[i] <- 1
       numdf1 <- (f(z+eps*e) - f(z))/eps
numdf1
}

numdf2 <- function(z,eps,i)  # 2nd order finite difference approximation to df at z
{
       e <- 0 * z
       e[i] <- 1
       numdf2 <- (f(z+eps*e) - f(z-eps*e))/(2.*eps)
numdf2
}


#
# as a check: this should give same answer as funcadm(x)
#
#
Rf <-  function(z){
  fortranf <- f(z)
  source(file.path(code,"readts.R"))
  err <- modts - obsts
  err[is.na(err)] <- 0
  out <- sum((err %*% t(w))^2)
  return(list(out=out,fortranf=fortranf,diff=out-fortranf))
}


#
# function to do latin hypercube sampling in hypercube bounder by [zlower:zupper]
#
dolhs <- function(nlhs, # number of function evaluations                 
                  zlo = zlower,
                  zhi = zupper
                 )
{
  dat <- randomLHS(nlhs,nz) # create sample of points in unit hypercube
  lhsz0 <- NA*dat           # initialise
  lhsz1 <- NA*dat           # initialise
  lhsf0 <- rep(NA,nlhs)
  lhsf1 <- rep(NA,nlhs)
  
  for (j in 1:nlhs)
  {
    lhsz0[j,] <- zlo + (zhi - zlo)*dat[j,] # scale starting points to fit in desired hypercube
  }

  for (ii in 1:nlhs)
  {
     print(paste("starting nlhs loop",ii))
     z <- lhsz0[ii,]
     source("findopt.R")
     print(ii)
     lhsz1[ii,] <- z1
     lhsf0[ii] <- f0
     lhsf1[ii] <- f1
     print("DONE in nlhs loop")
  }

  return(list(lhsz0 = lhsz0,
              lhsz1 = lhsz1,
              lhsf0 = lhsf0,
              lhsf1 = lhsf1))
}


#
# compare timeseries from two z vectors
#
zcompare <- function(z0=z0, # first z vector
                     z1=z1  # second z vector
                    )
{
  f(z0)
  source("readts.R")
  assign("modts0",modts,envir=.GlobalEnv)
  assign("err0",modts-obsts,envir=.GlobalEnv)

  f(z1)
  source("readts.R")
  assign("modts1",modts,envir=.GlobalEnv)
  assign("err1",modts-obsts,envir=.GlobalEnv)

}

rmse <- function(error=err1, # error time series
                 its=1       # observable
                )
{
out=sqrt(mean((error[,its])^2, na.rm=TRUE))
out
}


plotvar <- function(ivar=14:17,profile=1,type="p",imin=1,imax=nrow(R_out[[profile]]))
{

  ymin <- min(R_out[[profile]][imin:imax,ivar],na.rm=TRUE)
  ymax <- max(R_out[[profile]][imin:imax,ivar],na.rm=TRUE)

#  main=paste("Variable ",ivar,": ",names(R_out[[profile]])[ivar],sep="")

  plot(R_out[[profile]]$datetime[imin:imax],R_out[[profile]][imin:imax,ivar[1]],
       type=type,
       col=1,
       cex=0.3,
       ylim=c(ymin,ymax),
       xlab="Date",
       ylab=""
#       main=main
       )
  if (length(ivar) >= 2)
  {
     for (j in 2:length(ivar))
     {
     lines(R_out[[profile]]$datetime[imin:imax],R_out[[profile]][imin:imax,ivar[j]],
          type=type,
          col=j)
     }
  }
  legend("topright",
         legend=paste(ivar,": ",names(R_out[[profile]])[ivar],sep=""),
         col=1:length(ivar),
         pch=1)  
  axis(side=3,
       at=R_out[[profile]]$datetime[seq(from=imin,to=imax,length.out=8)], 
       labels=floor(seq(from=imin,to=imax,length.out=8)),
       tick=FALSE,
       line=-1
       )
  mtext("imin:imax",side=3,line=0.75) 
}


plotcompvar <- function(ivar=c(16,17),profile=1)
{

  plot(R_out[[profile]][,ivar[1]],R_out[[profile]][,ivar[2]],
       cex=0.5,
       xlab=names(R_out[[profile]])[ivar[1]],
       ylab=names(R_out[[profile]])[ivar[2]],
       col="black",           
       main=paste("Correlation between variables ",names(R_out[[profile]])[ivar[1]]," and ",names(R_out[[profile]])[ivar[2]],sep=""))
           

}








plotcorr <- function(its1=2,its2=3)
{
  imin <- min(c(err[,its1]), na.rm=TRUE)
  imax <- max(c(err[,its1]), na.rm=TRUE)

  ymin <- min(c(err[,its2]), na.rm=TRUE)
  ymax <- max(c(err[,its2]), na.rm=TRUE)

  plot(x = err[,its1],
       y = err[,its2],  
       ylim=c(ymin,ymax),
       xlim=c(imin,imax),
       type="p",
       xlab=paste("errors in timeseries #",its1,": ",tsnames[its1],sep=""),
       ylab=paste("errors in timeseries #",its2,": ",tsnames[its2],sep=""),
       main=paste("Error correlations for timeseries",tsnames[its1],"and",tsnames[its2]),
       col="red")
  abline(h=0,lty=2)
  abline(v=0,lty=2)     
}


















plottserr <- function(its=1,imin=1,imax=nrow(modts))
{

  ymin <- min(c(err[imin:imax,its]),na.rm=TRUE)
  ymax <- max(c(err[imin:imax,its]),na.rm=TRUE)

  plot(ts_time[imin:imax],err[imin:imax,its],
       type="b",
       cex=0.3,
#       xlim=c(imin,imax),
       ylim=c(ymin,ymax),
       xlab="Date",
       ylab=tsnames[its],
       main=paste("Errors in timeseries #",
                  its,
		  ": ",
		  tsnames[its],
		  sep=""),
       col="red")

  abline(h=0,lty=2)
  legend("topright",
         legend=c("model error"),
         col=c("red"),
         pch=c(19))

  axis(side=3,
       at=ts_time[seq(from=imin,to=imax,length.out=8)], 
       labels=floor(seq(from=imin,to=imax,length.out=8)),
       tick=FALSE,
       line=-1
       )
  mtext("imin:imax",side=3,line=0.75)

}








plotcompare <- function(its=2)
{
  aimin <- min(min(obsts[,its],na.rm=TRUE),na.rm=TRUE)
  aimax <- max(max(obsts[,its],na.rm=TRUE),na.rm=TRUE)

  plot(obsts[,its],
       modts[,its],
       cex=0.5,
       xlab=paste("Observed: ",tsnames[its]),
       ylab=paste("Modelled: ",tsnames[its]),
       col="black",
       xlim=c(aimin,
              aimax),
       ylim=c(aimin,
              aimax),              
       main=paste("Compare model and observations for timeseries #",its,": ",tsnames[its],sep=""))

rmse<-sqrt(mean(err[,its]^2,na.rm=TRUE))
corr<-cor.test(modts[,its],obsts[,its])$estimate
meanobs<-mean(obsts[,its],na.rm=TRUE)
meanmod<-mean(modts[,its],na.rm=TRUE)
fit1<-lm(modts[,its]~obsts[,its])

            
legend("bottomright",c(paste("RMSE: ",format(rmse,digits=4)),paste("R-squared:",format(corr,digits=4)),paste("obs mean: ",format(meanobs,digits=4)),paste("mod mean: ",format(meanmod,digits=4))))      
abline(0,1, lty=2, col="red")
abline(fit1$coefficients[1],fit1$coefficients[2], lty=2, col="blue")
legend("topleft",c("1:1 line","lm/best fit line"), col=c("red","blue"),lty=2)

}
if (!exists("runs"))
{
runs<-list()
}
storerun <-function(run_number=1)
{
runs[[run_number]]<<-list(obsts,modts,R_out)
}

chooserun <-function(run_number=1)
{
obsts<<-runs[[run_number]][[1]]
modts<<-runs[[run_number]][[2]]
R_out<<-runs[[run_number]][[3]]
}

whatvars<-function(profile=1)
{
print(paste("Available variables for profile ",profile,":"))
print(names(R_out[[profile]]))
}


comprunvar <- function(ivar=14:17,profile=1,run=1,imin=1,imax=nrow(R_out[[profile]]))
{

  ymin <- min(R_out[[profile]][imin:imax,ivar],na.rm=TRUE)
  ymax <- max(R_out[[profile]][imin:imax,ivar],na.rm=TRUE)

#  main=paste("Variable ",ivar,": ",names(R_out[[profile]])[ivar],sep="")

  plot(R_out[[profile]]$datetime[imin:imax],R_out[[profile]][imin:imax,ivar[1]],
       type="p",
       col=1,
       pch=1,
       cex=0.3,
       ylim=c(ymin,ymax),
       xlab="Date",
       ylab=""
#       main=main
       )
  lines(runs[[run]][[3]][[profile]]$datetime[imin:imax],runs[[run]][[3]][[profile]][imin:imax,ivar[1]],
          type="p",
          pch=2,
          col=1)
  
  if (length(ivar) >= 2)
  {
     for (j in 2:length(ivar))
     {
     lines(R_out[[profile]]$datetime[imin:imax],R_out[[profile]][imin:imax,ivar[j]],
          type="p",
          pch=1,
          col=j)
     lines(runs[[run]][[3]][[profile]]$datetime[imin:imax],runs[[run]][[3]][[profile]][imin:imax,ivar[j]],
          type="p",
          pch=2,
          col=j)
     }
  }
  legend("topright",
         legend=paste(ivar,": ",names(R_out[[profile]])[ivar],sep=""),
         pch=3,
         col=1:length(ivar))
  legend("topleft",
         legend=c("current run",paste("run",run)),
         pch=1:2)  
  axis(side=3,
       at=R_out[[profile]]$datetime[seq(from=imin,to=imax,length.out=8)], 
       labels=floor(seq(from=imin,to=imax,length.out=8)),
       tick=FALSE,
       line=-1
       )
  mtext("imin:imax",side=3,line=0.75) 
}









plotparamscatter <- function(choice="data"){
switch(choice,
       "ica"  = splom(S,       pch=".",main="Parameter scatter plot (ICA components)"),
       "data" = splom(zsample, pch=".",main="Parameter scatter plot"),
       "prc"  = splom(prc$x,   pch=".",main="Parameter scatter plot (Principal components)")
)
}
#
#

# now define some plotting functions
#

plotprc <- function(i1=1,
                    i2=2
                   ){
   biplot(prc, choices=c(i1,i2)) # a principal component plot
}

plotcorr <- function(){
corrplot(cor(zsample),method="ellipse",diag=FALSE,type="full")
}

plotcpdag <- function(){
   if (require(Rgraphviz)) {
       plot(pc.fit, main = paste("Estimated parameter CPDAG"))
   }
}

#
# define locf() - a local quadratic approximation to f() in the vicinity of optimum z1
#
locf <- function(z){
   locf <- f1 + as.numeric(t(df1) %*% (z-z1)) + 0.5 * as.numeric( t(z-z1) %*% hess %*% (z-z1) )
   locf
} 


gfs <- function(z){
   if (any(z < zlower) | any(z > zupper)){
       gfs <- sqrt(.Machine$double.xmax) } else {
       gfs <- sum((df(z)*(zupper-zlower))^2) }
       return(gfs)
}


# f1 and f2 could ne used with nlm()
#
f1 <- function(z){
   res <- f(z)
   attr(res,"gradient") <- df(z)
   res
}


f2 <- function(z){
   res <- f(z)
   attr(res,"gradient") <-  df(z)
   attr(res,"hessian")  <- d2f(z)
   res
}

#
# define preconditioned versions of f and df
#
# we have y = S %*% z
#
# with S = sqrt|D| %*% t(L)
#
# where hess  = L %*% D %*% t(L) from gchol()
g <- function(y){
     return(f(Sinv%*%y))
}

dg <- function(y){
return(Sinv %*% df(Sinv %*% y))
}


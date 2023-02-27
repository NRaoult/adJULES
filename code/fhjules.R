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
# R code to initialise the adjoint model of JULES, with Hessian.
# There is a option to change to the vector tangent linear JULES version
#
# (mabye rename as load_adjules.R?)

source(file.path(code,"setlibraries.R"))
#
# define objective function
#
f <- function(z,xoff)
{
   x <- xoff       # default value
   x[xvary] <- z
   f <- func(x)
   f
}

df <- function(z,xoff)
{
   x <- xoff       # default value
   x[xvary] <- z
   df <- dfunc(x)[xvary]
   names(df) <- znames
   df
}
# define hessian function - only when using fhjules  

d2f <- function(z,xoff=xoff){
   x <- xoff       # default value
   x[xvary] <- z
   f <- d2func(x)[xvary,xvary]
   f
}   

# check this for global and local variables
setupFort <- function(){
          .Fortran("setwindow",Rwindow=window)    
          .Fortran("setw",w=w)
          .Fortran("setsmooth",lsmooth=lsmooth)    # choose logical value for smoothing
          .Fortran("setmask",lmask=lmask)          # choose whether to mask nighttime fluxes
          .Fortran("setPM",M=M)
          .Fortran("setxvary",xvaryin=xvary)
}


func  <- function(x, type ='adm'){  #other tl
          
          setupFort()
          if ( (sum(x>xmax) + sum(x < xmin)) > 0 ){ print("x outside bounds in funcadm - resetting") }

          xlocal <- pmin(pmax(x,xmin),xmax)

          temp <- .Fortran(paste0("vfunc_wrapper_",type),x1=xlocal,v=double(m))        
          func <- as.double(temp$v) 
    if ( exists("rec",envir=.GlobalEnv)){assign("rec",rbind(rec,c(xlocal,func)),envir=.GlobalEnv)}
    if (!exists("rec",envir=.GlobalEnv)){assign("rec",c(xlocal,func),envir=.GlobalEnv)}
    func
}

#
# define gradient function
#	  
dfunc <- function(x){#,type = 'adm'){   #vtl
          
          setupFort()
          if ( (sum(x>xmax) + sum(x < xmin)) > 0 ){ print("x outside bounds in func - resetting") }

          xlocal <- pmin(pmax(x,xmin),xmax)

          temp <- .Fortran("dfuncadm",x=xlocal,grady=double(n))
          dfunc <- as.double(temp$grady) 
	  dfunc
	  }


d2func <- function(x){
          setupFort()
          temp <- .Fortran("d2func_wrapper",x=x,hess=matrix(0,nrow=n,ncol=n))
          d2func <- matrix(as.double(temp$hess),ncol=n,nrow=n) 
          d2func
} 
#
# define function to perform z1 $ opt$parisation
#
minim <- function(z,hessian=TRUE,parscale=zupper-zlower){
    check1 <- sum(z < zlower)
    check2 <- sum(z > zupper)
    if (check1+check2 > 0) print("minim WARNING: z is not within [lower:upper]")
  
    minim <- optim(par     = z,   
                   fn      = f,
                   gr      = df, 
                   method  = "L-BFGS-B",
                   lower   = zlower,
                   upper   = zupper, 
                   hessian = hessian,
                   control = list(REPORT=1,
                                  trace=1,          # trace=1 usually, trace=6 very complete
                                  parscale=parscale
                                  )
                   )
    minim
}

#dyn.load("vtjules.so") 
#dyn.load("adjules.so")  
dyn.load(file.path(code,"fhjules.so"))  

temp <- .Fortran("setfunc",n=integer(1),m=integer(1))
  n <- temp$n
  m <- temp$m
rm(temp)

M <- mat.or.vec(n,n)   # initialise background term

source(file.path(code,"define_functions.R"))

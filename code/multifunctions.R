# loads in multisite_functions
multif <- function(zlocal,lgetw=FALSE){
   sum <- 0
   for (ilocal in 1:nsites){

      unlink(  c("modts.da" , "obsts.da" , "control.jin", "avmodts.da", "avobsts.da")  )                   
      temp1 = file.copy( multidatfile[ilocal] , "obsts.da" ) ; if (!temp1) stop("Error: observations not copied over")
      temp2 = file.copy( multijinfile[ilocal],"control.jin" ); if (!temp2) stop("Error: jinfile not copied over")
      rm(temp1 , temp2)

      .Fortran("initfunc",n=n,x=double(n))
      if (lgetw){f(zlocal,xoffmulti[,ilocal]); source(file.path(code,"getw.R"));multiw[,,ilocal] <<- w}
      assign("w",multiw[,,ilocal],envir=.GlobalEnv)
      sum <- sum + f(zlocal,xoffmulti[,ilocal]) 

      .Fortran("postfunc", n=n, x=xoffmulti[,ilocal],  m=m, y=double(m))
   }
   sum
}



multidf <- function(zlocal){
   sum <- 0
   for (ilocal in 1:nsites){
      unlink(  c("modts.da" , "obsts.da" , "control.jin", "avmodts.da", "avobsts.da")  )                   
      temp1 = file.copy( multidatfile[ilocal] , "obsts.da" ) ; if (!temp1) stop("Error type CARROT in setup.R")
      temp2 = file.copy( multijinfile[ilocal],"control.jin" ); if (!temp2) stop("Error type SPUD in setup.R")
      rm(temp1 , temp2)
      .Fortran("initfunc",n=n,x=double(n)) #          .Fortran("postfunc",n=n,m=m,x=x,y=double(m))
      assign("w",multiw[,,ilocal],envir=.GlobalEnv)
      sum <- sum + df(zlocal,xoffmulti[,ilocal])  
      .Fortran("postfunc", n=n, x=xoffmulti[,ilocal],  m=m, y=double(m))
   }
   sum
}



multid2f <- function(zlocal){
   sum <- 0
   for (ilocal in 1:nsites){
      unlink(  c("modts.da" , "obsts.da" , "control.jin", "avmodts.da", "avobsts.da")  )                   
      temp1 = file.copy( multidatfile[ilocal] , "obsts.da" ) ; if (!temp1) stop("Error type CARROT in setup.R")
      temp2 = file.copy( multijinfile[ilocal],"control.jin" ); if (!temp2) stop("Error type SPUD in setup.R")
      rm(temp1 , temp2)
 
      .Fortran("initfunc",n=n,x=double(n)) #
      assign("w",multiw[,,ilocal],envir=.GlobalEnv)
      sum <- sum + d2f(zlocal,xoffmulti[,ilocal])  

      .Fortran("postfunc", n=n, x=xoffmulti[,ilocal],  m=m, y=double(m))
   }
   sum
}


#
# define function to perform minimisation
#
multiminim <- function(z){
     check1 <- sum(z < zlower) ; check2 <- sum(z > zupper)
     if (check1+check2 > 0) print("warning: z is not within [lower:upper]")
  
     z <- pmin(pmax(z,zlower),zupper) # Jupp added to prevent problems
     parscale <- zupper-zlower

     minim <- optim(par     = z,
                    fn      = multif,
                    gr      = multidf, 
                    method  = "L-BFGS-B",
                    lower   = zlower,
                    upper   = zupper, 
                    hessian = FALSE,
                    control = list(REPORT= 1,
                                   trace = 1,         # trace=1 usually, trace=6 very complete
                                   maxit = 150,       # can change this 
                                   factr = 1.0e9,     # Luke changed to eliminate convergence problem
                                 parscale= parscale
                                  )
                   )
     minim
}

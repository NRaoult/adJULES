# *****************************COPYRIGHT*******************************
# (c) University of Exeter 2017 . All rights reserved.
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
# R code to calculate the weight matrix w
#
#f(z)         # get cost with no parameter cost to save timeseries to disk

source(file.path(code,"readts.R"))
err <- modts - obsts
err[is.na(err)] <- 0
npoints <- length(which(err[,tsuse]!=0)) # number of time points contributing to cost. Replaces nt*ntsuse.

# weight background term
print(paste("setting M with alpha",alpha))
M <- alpha * sqrt(npoints-(nz/nsites)) * diag(1./(xmax-xmin))

#avoid singular matrices
temp        <- t(err) %*% err
mask        <- diag(temp)==0
tsuse[mask] <- FALSE

ntsuse <- sum(tsuse)

c1 <- ntsuse*(t(err[,tsuse]) %*% err[,tsuse])  # non central 2nd moment of errors

w <- matrix(0,nts,nts) # default

w[tsuse,tsuse] <- chol(solve(c1))
 
print(paste("Data points to fit:",npoints))
print(paste("Parameters to vary:",nz))
print(paste("Degrees of freedom:",npoints-nz)) # need to change this for multisite

if (npoints > nz) {
   w <- sqrt(npoints-(nz/nsites)) * w  # estimate covariance for chi-square modelling
} else {
   print("WARNING: more parameters than data points!!!")
   print("consider shorter averaging window in order to make more data points")
}

colnames(w) <- tsnames
rownames(w) <- tsnames

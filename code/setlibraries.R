
# *****************************COPYRIGHT*******************************
# (c) University of Exeter 2014 . All rights reserved.
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
# script to install required R libraries
#

# may need to install non-cran libraries as follows@
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("Rgraphviz")
#   biocLite("graph")
#   biocLite("RBGL")

library(corrplot)
#library(pcalg)
library(tmvtnorm) 
library(ellipse)
#library(bdsmatrix)    # required for calls to gchol()
#library(lattice)
#library(lhs)
library(corpcor)      # required for calls to make.positive.definite()
library(KernSmooth)

library('R.utils')

#library(optimx)       # loads many optimisation routines for use in optimx()
#library(emulator)
#library(trustOptim)   # required for calls to trust.optim()

# library(minqa)      # required for calls to bobyqa()
library(RColorBrewer)
#library(colorspace)



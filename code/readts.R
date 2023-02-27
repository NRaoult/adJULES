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
# R code to read in binary timeseries data for plotting. Binary timeseries files genereated by running
# f(z) and changing the 'window' before running will give different averaging
#
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

fullmodts <- read_timeseries('modts.da') # full modelled timeseries
fullobsts <- read_timeseries('obsts.da') # full observed timeseries
modts <- read_timeseries('avmodts.da') # averaged modelled timeseries
obsts <- read_timeseries('avobsts.da') # averaged observed timeseries

#!/usr/bin/env Rscript

#args <- commandArgs(trailingOnly=TRUE)
#print(args)

source('../../set_paths.R')

##########################
#
#  EDIT HERE to define experiment
#
#########################
#
# When defining sets, a few tips:
# - 1:nrow(data) will multisite over the whole set
# - sample(1:nrow(data),5) will multisite over a randomly selected 5 sites
# - for NT, the set can be split boreal=c(1:11,17,25,26), rest=c(12:16,18:24,27:35)
# - for BT, decidious sites are 1:18, evergreen are 19:28

experiment   <- 'EXPERIMENT'
type	     <- 'SUBSET' 
alpha        <- 'ALPHA'
update_start <- 'ON'
validation   <- 'off'

stream     <- 'GPP'
window     <- 17520/12 #180.0
windtype   <- 'monthly'# 'daily'

source(file.path(code,'load_experiment_data.R'))

# averaged out z (for clustering at the same start)
z_av <- c(0.0544,  0.0840,  0.8700,  0.6000, 35.8000, 1.1000,  0.05, 0.0850) #canht = 7.688, dlai =  0.0500

ss<-''
# set up the different 'sets' to simultaneously optimise over
if (type == 'all'){
	set <- 1:nrow(data)
} else if (type == 'boreal'){
	set <- c(1:11,17,25,26)
} else if (type == 'rest'){
	set <- c(12:16,18:24,27:35)
} else if (type == 'deci'){
	set <- c(1:18)
} else if (type == 'ever'){
	set <- c(19:28)
} else if (type == 'ran5'){
	set <- sample(1:nrow(data),5)
        ss <- toString(set)
} else if (type == 'ran5deci'){
        set <- sample(1:18,5)
	ss <- toString(set)
} else if (type == 'single'){     #running single site through multisite code (same)
	for (ss in 1:nrow(data)){
		set <- ss
		if (alpha == 'find'){
 	            source(file.path(adjules,'alpha_loop.R'))
 	        } else if (alpha == 'brute'){
 	            for (alpha in c(0,1))  source(file.path(code,'multisite.R')) #seq(0,1,0.1)
		} else if (alpha == 'sqrt'){
			alpha <- 1/sqrt(2)
			source(file.path(code,'multisite.R'))
 	        } else {
                    alpha <- as.numeric(alpha)
                    source(file.path(code,'multisite.R'))
                }
	}
} else {
	set <- scan(text=type,sep=',',what =numeric(),quiet=TRUE)
	type <- 'chosen'
	ss <- toString(set)
}

if (type != 'single'){#ss <-''
	if (alpha == 'find'){
 	            source(file.path(adjules,'alpha_loop.R'))
 	        } else if (alpha == 'brute'){
 	            for (alpha in c(0,1))  source(file.path(code,'multisite.R')) #seq(0,1,0.1)
			} else if (alpha == 'sqrt'){
			alpha <- 1/sqrt(2)
			source(file.path(code,'multisite.R'))
 	        } else {
                    alpha <- as.numeric(alpha)
                    source(file.path(code,'multisite.R'))
                }
}

# Script to load the data in experiment data, which contains information about all the sites needed to load in jinfile and data data files
# 

### switches for different experiments
# Change experiment year to validation year (might be better way to switch this)
stream <- 'GPP'
verifydir = file.path(code,'../jules_verify_data',experiment)

### Load experiment data 
data <- as.matrix(read.table(paste0(verifydir,'/experiment_data_',experiment,'.txt'), header=TRUE))

locs  <- data[,1]
start <- data[,2]
PFT   <- data[,4]
exp_year <- data[,5]

PFT_type <- c('BT','NT','C3','C4','Sh')
nPFT     <- length(PFT_type)
pft      <- which(PFT_type %in% experiment)

nts <- 6                                              # the number of timeseries variables
tsnames <- c("NEE","H","LE","Tstar","GPP","Resp")     # names of timeseries variables
	if (stream == 'GPP') tsuse <- c(FALSE,FALSE, TRUE, FALSE, TRUE,  TRUE)  
	if (stream == 'NEE') tsuse <- c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) 
names(tsuse) <- tsnames
w <- diag(as.integer(tsuse))                          # a dummy value to get started 

lmask   <- 2
lsmooth <- 1

## Load in parameter data
param_data  <- read.table(file.path(code,'parameters.txt'), sep = '', header=TRUE) #
xnames      <- param_data[,1]
xmin        <- param_data[,2]
xmax        <- param_data[,3]
names(xmax) <- xnames
names(xmin) <- xnames

parameterchoice <- c(2,7,12,17,22,46,68,73) + pft -1            # choose common parameters ,68
# nl0, alpha,f0, tlow, tupp, rootd_ft, dcatch_dlai (68), dqcrit (canht =32)

xvary <- rep(FALSE,length(xnames))
names(xvary) <- xnames
xvary[parameterchoice] <-TRUE 

znames        <- xnames[xvary]
zlower        <- xmin[xvary]
zupper        <- xmax[xvary]
names(zlower) <- znames
names(zupper) <- znames

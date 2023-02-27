#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
# EDIT here if different path
pwd <- getwd()

# Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}
 
# Help section - need to figure out when alpha defined but not subset
if("--help" %in% args) {
  cat("
      The create temp folder for Script
 
      Arguments:
      arg1   - name of temp folder            (NO default) 
      arg2   - PFT ('BT','NT','C3','C4','Sh') (NO default)
      arg3   - subset of sites (all,single, boreal,ever,deci,rest,ran5) (default=all)
      arg4   - alpha, either input a value or search for optimal with 'find' or 'brute' (default=find)
      arg5   - update starting point (default=off)
      help   - print this text
 
      Example:
      ./create_temp.R deci_loop BT deci \n\n")
 
  q(save="no")
}

if(args[3]=='set'){
readline("Input set of sites") -> ichoice
args[3] <- ichoice
}

## Arg3 default
if(is.na(args[3])) {
  args[3] <- 'all'
}
 
## Arg4 default
if(is.na(args[4])) {
  args[4] <- 'find'
}

## Arg5 default
if(is.na(args[5])) {
  args[5] <- 'off'
}

exname     <- args[1]
experiment <- args[2]
type       <- args[3]
alpha      <- args[4]
update     <- args[5]

dir_name <- paste0('temp/temp_',exname)
system(paste0('mkdir -p ',dir_name,'/{OUTPUT,jules_verify_data}'))
system(paste0('cp -r ',pwd,'/jules_verify_data/ temp/temp_',exname,'')) #need to change paths in jin files before removing this line

values <- matrix(c( 'EXPERIMENT', experiment,
        	          'SUBSET',     type,
           	        'ALPHA',      alpha,
		                'ON',         update, 
                  ), ncol=2,byrow=TRUE)

x <- readLines(paste0(pwd,'/run_template.R'))

for (i in 1:nrow(values)) x <- gsub(values[i,1],values[i,2], x)
cat(x,file=paste0(pwd,'/temp/temp_',exname,'/run.R'),sep='\n',append=FALSE)

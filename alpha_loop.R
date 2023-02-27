#
# A script to find alpha using a bisection algorithm. Slow but garanteed to converge compared to other method.
# Expect 7 iterations (log2(1/0.01)), though premature rounding might extend this to a couple more.
#
# Assumption made: linearity. Clear step between positive definite and not positive definite
#
# Other options to consider: - Secant method doesnt garantee convergence, though regula falsi method could
#                            - Using eigenvalues to define function


#file_name <- paste0(result,'alpharange_',experiment,'_',type,'.txt')
file_name <- paste0('alpharange_',experiment,'_',type,'_',window,'.txt')
write(paste('Experiment',experiment,'using',type,'as subset, has alpha in range'), file_name)

# A script to run through different alpha values, narrowing the range for positive definite

low_bound <- 0
high_bound <- 1
maxcount <- 10

alpha <- 0
source(file.path(code,'multisite.R'))
alpha <- 1
source(file.path(code,'multisite.R'))

while (high_bound-low_bound > 0.02){ # & i < maxcount){
	#i = i + 1
	if (is.positive.definite(hess_taf)){
		high_bound <- alpha		
		alpha <- round((high_bound-low_bound)/2+low_bound,2)	
		if (alpha <0 ) break
		source(file.path(code,'multisite.R'))
	} else {
		if (alpha == high_bound) {high_bound <- round(alpha*1.5,2)}
		low_bound <- alpha
		alpha <- round((high_bound-low_bound)/2+low_bound,2)
		if (alpha <0 ) break
		source(file.path(code,'multisite.R'))
	} 
	write(paste(low_bound,':',high_bound),file_name,append=TRUE)
}

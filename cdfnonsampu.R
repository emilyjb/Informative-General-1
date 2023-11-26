

cdfnonsampuOld <- function(u, lambda0, lambda1, sigma2u, tau, r){
	esw <- exp(lambda0 + tau^2/2 + sigma2u*lambda1^2/2)
	num <- esw*pnorm(u, mean = sigma2u*lambda1, sd = sqrt(sigma2u)) - pnorm(u, mean = 0, sd = sqrt(sigma2u))
	den <- esw - 1
	num/den -r
}

genudnonsamp <- function(lambda0, lambda1, sigma2u, tau){ 
	r <- runif(1, 0, 1);  
	uniroot(cdfnonsampu, c(-10000, 10000), lambda0 =  lambda0 ,lambda1 = lambda1, sigma2u = sigma2u, tau = tau, r=r)$root 
}

cdfnonsampu <- function(u, lambda0, lambda1, sigma2u, tau, r){
	term1  <-  exp(lambda0 + tau^2/2 + sigma2u*lambda1^2/2 - log(exp(lambda0 + tau^2/2 + sigma2u*lambda1^2/2)-1))
	term2 <- term1 - 1
	num <- term1*pnorm(u, mean = sigma2u*lambda1, sd = sqrt(sigma2u)) - term2*pnorm(u, mean = 0, sd = sqrt(sigma2u))
	num  -r
}



transfun <- function(x){
	 qnorm(pchisq(x, df = 1))
}

invtransfun <- function(y){
	qchisq(pnorm(y), df = 1)
}

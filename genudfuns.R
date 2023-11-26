

compnumlogist <- function(u, beta, xs, ys, areafacsamp, d){
	xsd <- xs[areafacsamp == d]
	ysd <- ys[areafacsamp == d]
	etad <- as.vector(cbind(1, xsd)%*%beta + u)
	pd <- exp(etad)/(1 + exp(etad))
	a1 <- prod( (pd^ysd)*( (1 - pd)^(1-ysd)))
      a2 <-  exp( sum(  (log(pd)*ysd)+ ( log((1 - pd))*(1-ysd)) +100 ))
	ifelse(!is.na(a1), a1, a2)

}


genudlogist <- function(d, sigma2u, beta, xs, ys, areafacsamp, L = 200){
	uL <- rnorm(L, mean = 0, sd = sqrt(sigma2u ))
	probnumd <- sapply(uL, compnumlogist, beta, xs, ys, areafacsamp, d)
	upickd <- sample(uL, size= 1, prob=probnumd/sum(probnumd))
	upickd
}



compnumlogistMod <- function(u, beta, xs, ys, areafacsamp, sigma2u, uhatmean, uhatvar, d){
        xsd <- xs[areafacsamp == d]
        ysd <- ys[areafacsamp == d]
        etad <- as.vector(cbind(1, xsd)%*%beta + u)
        pd <- exp(etad)/(1 + exp(etad))
        prod( (pd^ysd)*( (1 - pd)^(1-ysd)))*dnorm(u, mean = 0, sd = sqrt(sigma2u))/dnorm(u, mean = uhatmean[d], sd = sqrt(uhatvar[d]))
}

genudlogistMod <- function(d, sigma2u, beta, xs, ys, areafacsamp, uhatvar, uhatmean, L = 200){
        uL <- rnorm(L, mean = 0, sd = sqrt(uhatvar[d]))
        probnumd <- sapply(uL, compnumlogistMod, beta, xs, ys, areafacsamp, sigma2u, uhatmean, uhatvar, d)
        upickd <- sample(uL, size= 1, prob=probnumd/sum(probnumd))
        upickd
}




compufunareapar <- function(uvec, ws, ys, xs, beta, sigma2u, areafacsamp, areapar){
	names(uvec) <- sort(unique(areafacsamp))
	eta <-as.vector( cbind(1, xs)%*%beta + uvec[as.character(areafacsamp)])
	pij <- exp(eta)/(1 + exp(eta))
	prodpij <- tapply( (pij^ys)*( (1 - pij)^(1 - ys)), areafacsamp, prod)
	1/sqrt(areapar[3])*dnorm( ( log(ws) - areapar[1] - areapar[2]*uvec)/sqrt(areapar[3]) )*prodpij
}




compufunareapard <- function(u , ws, ys, xs, beta, sigma2u, areafacsamp, areapar, d ){
	ysd <- ys[areafacsamp == d]
	xsd <- xs[areafacsamp == d]
	eta <-as.vector( cbind(1, xsd)%*%beta + u)
	pij <- exp(eta)/(1 + exp(eta))
	prodpij <- prod ( (pij^ysd)*( (1 - pij)^(1 - ysd)) )
	d2 <- which( sort(unique(areafacsamp)) == d)
	1/sqrt(areapar[3])*dnorm( ( log(ws[d2]) - areapar[1] - areapar[2]*u )/sqrt(areapar[3]) )*prodpij
}




objfunareapars <- function(areapar, ws, ys, xs, beta, sigma2u, areafacsamp, T = 200){
	ubigmat <- matrix( rnorm(length(unique(areafacsamp))*T, mean = 0, sd = sqrt(sigma2u)), nrow = length(unique(areafacsamp)), ncol = T)
	compuout <- apply(ubigmat, 2, compufunareapar,  ws, ys, xs, beta, sigma2u, areafacsamp, areapar)
	-sum(log(apply(compuout, 1, mean)	))
}



objfunareaparsGH <- function(areapar, ws, ys, xs, beta, sigma2u, areafacsamp, T = 20){
	ghp <- gauss.quad.prob(T, dist = "normal", mu = 0, sigma  = sqrt(sigma2u))
	wts <- ghp$weights
	nds <- ghp$nodes
	ubigmat <- kronecker( rep(1,length(unique(areafacsamp))), t(nds))
	wbigmat <- kronecker( rep(1,length(unique(areafacsamp))), t(wts))
	compuout <- apply(ubigmat, 2, compufunareapar,  ws, ys, xs, beta, sigma2u, areafacsamp, areapar)
	-sum(log(apply(compuout*wbigmat, 1, sum)	))
}


 
objfunareaparsGHbprime <- function(areapar, ws, ys, xs, beta, sigma2u, areafacsamp, T = 20){
	ghp <- gauss.quad.prob(T, dist = "normal", mu = 0, sigma  = sqrt(sigma2u))
	wts <- ghp$weights
	nds <- ghp$nodes
	ubigmat <- kronecker( rep(1,length(unique(areafacsamp))), t(nds))
	wbigmat <- kronecker( rep(1,length(unique(areafacsamp))), t(wts))
	compuout <- apply(ubigmat, 2, compufunareaparbprime ,  ws, ys, xs, beta, sigma2u, areafacsamp, areapar)
	-sum(log(apply(compuout*wbigmat, 1, sum)	))
}


compufunareaparbprime <- function(uvec, ws, ys, xs, beta, sigma2u, areafacsamp, areapar){
	names(uvec) <- sort(unique(areafacsamp))
	eta <-as.vector( cbind(1, xs)%*%beta + uvec[as.character(areafacsamp)])
	pij <- exp(eta)/(1 + exp(eta))
	prodpij <- tapply( (pij^ys)*( (1 - pij)^(1 - ys)), areafacsamp, prod)
	muw <- exp(areapar[1] + areapar[2]*uvec)/(1 + exp(areapar[1] + areapar[2]*uvec))
	norm <- 1/beta( (1 - muw)*areapar[3], muw*areapar[3] + 1)
	 norm*(ws - 1)^( (1 - muw)*areapar[3] -1)*ws^(-areapar[3]-1)*prodpij
}


 



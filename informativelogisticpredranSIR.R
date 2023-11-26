
betahat <- fixef(glmersamp)

r <- 0
R <- 100

ybarinfrs <- c()
ybarnoninfrs <- c()
yvarinfrs <- c()
yvarnoninfrs <- c()
yoddsinfrs <- c()
yoddsnoninfrs <-  c()


repeat{

	r <- r + 1
	# For new output
	#ugenr <-   sapply(1:D, genudlogistgrid , sigma2u, betahat, xs, ys, areafacsamp, 100) 
	ugenr <-   sapply(1:D, genudlogist  , sigma2u, betahat, xs, ys, areafacsamp, 100) 

	names(ugenr) <- 1:D
	ugenlong <- ugenr[as.character(areafacpop)]
	# test
	#ugenlong <- ranefsamp[as.character(areafacpop)]
	#
	etapopInf <- as.vector(cbind(1,xN)%*%coefsamp + ugenlong + bhat)
	etapopNonInf <- as.vector(cbind(1,xN)%*%coefsamp + ugenlong )
	phatinf <- exp(etapopInf)/(1 + exp(etapopInf))
	phatnoninf <- exp(etapopNonInf)/(1 + exp(etapopNonInf))
	ygeninf <- sapply(phatinf, function(x){ rbinom(1, prob = x, size = 1)})
	ygennoninf <- sapply(phatnoninf, function(x){ rbinom(1, prob = x, size = 1)})

	ybarinfr <- tapply(ygeninf, areafacpop, mean)
	ybarnoninfr <- tapply(ygennoninf, areafacpop, mean)

	yvarinfr <- tapply(ygeninf, areafacpop, var)
	yvarnoninfr <- tapply(ygennoninf, areafacpop, var)

	#yoddsinfr <- log(tapply(1 - ygeninf, areafacpop, mean)/tapply( ygeninf, areafacpop, mean))
	#yoddsnoninfr <-log( tapply(1- ygennoninf, areafacpop, mean)/tapply( ygennoninf, areafacpop, mean))

	yoddsinfr <-  tapply(1 - ygeninf, areafacpop, mean)/tapply( ygeninf, areafacpop, mean) 
	yoddsnoninfr <-  tapply(1- ygennoninf, areafacpop, mean)/tapply( ygennoninf, areafacpop, mean) 
	
	ybarinfrs <- rbind(ybarinfrs, ybarinfr)
	ybarnoninfrs <- rbind(ybarnoninfrs, ybarnoninfr)
	yvarinfrs <- rbind(yvarinfrs, yvarinfr)
	yvarnoninfrs <- rbind(yvarnoninfrs, yvarnoninfr)
	yoddsinfrs <- rbind(yoddsinfrs, yoddsinfr)
	yoddsnoninfrs <- rbind(yoddsnoninfrs, yoddsnoninfr)

	if(r == R){break}

}

ybarinf <- apply(ybarinfrs, 2, mean)
ybarnoninf <- apply(ybarnoninfrs, 2, mean)
g1infbar <- apply(ybarinfrs, 2, var)

yvarinf <- apply(yvarinfrs, 2, mean)
yvarnoninf <- apply(yvarnoninfrs, 2, mean)
g1infvar <- apply(yvarinfrs, 2, var)

yoddsinf <- apply(yoddsinfrs, 2, mean)
yoddsnoninf <- apply(yoddsnoninfrs, 2, mean)
g1infodds <- apply(yoddsinfrs, 2, var)


ybarinfs <- rbind(ybarinfs, ybarinf)
ybarnoninfs <- rbind(ybarnoninfs, ybarnoninf)
g1infbars <- rbind(g1infbars, g1infbar)

yvarinfs <- rbind(yvarinfs, yvarinf)
yvarnoninfs <- rbind(yvarnoninfs, yvarnoninf)
g1infvars <- rbind(g1infvars, g1infvar)

yoddsinfs <- rbind(yoddsinfs, yoddsinf)
yoddsnoninfs <- rbind(yoddsnoninfs, yoddsnoninf)
g1infoddss <- rbind(g1infoddss, g1infodds)









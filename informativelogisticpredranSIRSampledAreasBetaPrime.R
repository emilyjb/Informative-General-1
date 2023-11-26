
betahat <- fixef(glmersamp)

r <- 0
R <- 100

ybarinfrs <- c()
ybarnoninfrs <- c()
yvarinfrs <- c()
yvarnoninfrs <- c()
yoddsinfrs <- c()
yoddsnoninfrs <-  c()

xNareasamp <- xN#[which(IDlong == 1)]

###### First do sapled areas:
repeat{

	r <- r + 1
	ugenr <- sapply(sort(unique(areafacsamp)), genudlogist , sigma2u, betahat, xs, ys, areafacsamp, 100)
	names(ugenr) <- sort(unique(areafacsamp))
	ugenlong <- ugenr[as.character(areafacpop[which(areafacpop %in% areafacsamp)])]
	etapopInf <- as.vector(cbind(1,xNareasamp)%*%coefsamp + ugenlong - bhat)
	etapopNonInf <- as.vector(cbind(1,xNareasamp)%*%coefsamp + ugenlong )
	phatinf <- exp(etapopInf)/(1 + exp(etapopInf))
	phatnoninf <- exp(etapopNonInf)/(1 + exp(etapopNonInf))

	ygeninf <- sapply(phatinf, function(x){ rbinom(1, prob = x, size = 1)})
	ygeninf[sampindex] <- ys
	ygennoninf <- sapply(phatnoninf, function(x){ rbinom(1, prob = x, size = 1)})
	ygennoninf[sampindex] <- ys

	ybarinfr <- tapply(ygeninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)
	ybarnoninfr <- tapply(ygennoninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)

	yvarinfr <- tapply(ygeninf, areafacpop[which(areafacpop %in% areafacsamp)], var)
	yvarnoninfr <- tapply(ygennoninf, areafacpop[which(areafacpop %in% areafacsamp)], var)

	#yoddsinfr <- log(tapply(ygeninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)/tapply(1-ygeninf, areafacpop[which(areafacpop %in% areafacsamp)], mean))
	#yoddsnoninfr <-log( tapply(ygennoninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)/tapply(1-ygennoninf, areafacpop[which(areafacpop %in% areafacsamp)], mean))

	yoddsinfr <- tapply(1-ygeninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)/tapply( ygeninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)
	yoddsnoninfr <-  tapply(1-ygennoninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)/tapply( ygennoninf, areafacpop[which(areafacpop %in% areafacsamp)], mean)
	
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

#checkoddsInf <- sum(yoddsinf == Inf) + sum(yoddsnoninf==Inf) + sum(g1infodds == Inf) > 0

#if(checkoddsInf){break}







B <- 100
b <- 0

ybarinfbs <- c()
g1infbarbs <- c()

yvarinfbs <- c()
g1infvarbs <-c()

yoddsinfbs <- c()
g1infoddsbs <-c()

ybarpopbs <- c()
yvarpopbs <-  c()
yoddspopbs <-  c()

parbs <- c()

repeat{

	b <- b + 1
	
	ub <- rnorm(D, mean = 0, sd = sqrt(sigma2u))
	names(ub) <- 1:D

	psb <- exp(coefsamp[1] + coefsamp[2]*xs + ub[as.character(areafacsamp)])/(1 + exp(coefsamp[1] + coefsamp[2]*xs + ub[as.character(areafacsamp)]))
	ysb <- rbinom(length(psb), size=1, prob = psb)

	#pcompInfN <- exp(coefsamp[1] + coefsamp[2]*xN  - bhat + ub[as.character(areafacpop)])/(1 + exp(coefsamp[1] + coefsamp[2]*xN  - bhat + ub[as.character(areafacpop)]))
	 
	#psampInfN <- exp( - bhat+ coefsamp[1] + coefsamp[2]*xN   + uball[as.character(areafacpop)])/(1 + exp(coefsamp[1] + coefsamp[2]*xN  + uball[as.character(areafacpop)]))
	#pcompInfN <- exp(coefsamp[1] + coefsamp[2]*xN  - bhat + uball[as.character(areafacpop)])/(1 + exp(coefsamp[1] + coefsamp[2]*xN  - bhat + uball[as.character(areafacpop)]))
	#denN <- 1 + exp(-areapar[1])*((1 - psampInfN) + psampInfN*exp(-areapar[2]))
	#phatinfN <- (psampInfN/denN + exp(-areapar[1])*((1 - psampInfN) + psampInfN*exp(-areapar[2]))*pcompInfN/denN	)
	#ypopNb <- rbinom(length(phatinfN), size=1, prob = phatinfN)

	pcb <- exp(coefsamp[1] + coefsamp[2]*xN - bhat + ub[as.character(areafacpop)])/(1 + exp(coefsamp[1] + coefsamp[2]*xN - bhat + ub[as.character(areafacpop)]))
	yNb <- rbinom(length(pcb), size=1, prob = pcb)
	yNb[sampindex] <- ysb
	
	#yNb[! which(areafacpop  %in% areafacsamp)] <- ypopNb[! which(areafacpop  %in% areafacsamp)]

	mupijb <- exp( ahat + bhat*ysb)/( 1 + exp(ahat + bhat*ysb))
	pijsb <- rbeta(length(mupijb), mupijb*phiijhat + 1, (1 - mupijb)*phiijhat)
	wijsb <- 1/pijsb
	#mupisb <- exp(areapar[1] + areapar[2]*ub)/(1 + exp(areapar[1] + areapar[2]*ub))
	#wisb <- 1/rbeta(length(mupisb), shape1 =  mupisb*areapar[3] + 1,  shape2 = (1 - mupisb)*areapar[3])	
	

	#wisb <- replicate(100000, 1/rbeta(length(mupisb), shape1 =  mupisb*areapar[3] + 1,  shape2 = (1 - mupisb)*areapar[3]) )
	#mh <- apply(wisb, 1, mean)
	#mt <- 1 + exp( - areapar[1] - areapar[2]*ub)

	#plot(mt, mh)

	glmersampb <- glmer(cbind(ysb, 1 - ysb) ~ xs + (1 | areafacsamp), family = binomial(link = "logit"), nAGQ = 20)
	lmwijb  <- lm(log(wijsb -1 )~  ysb )
	bhatb <- lmwijb$coef[2]
	ahatb <- lmwijb$coef[1]
	mupijinitb <- exp(-ahatb - bhatb*ysb)/(1 + exp(-ahatb - bhatb*ysb))
	phijinitb <- 1/mean((pijsb  - mupijinitb)^2/(mupijinitb*(1 - mupijinitb)))

	pijparb <- optim( c(-ahatb, -bhatb, phijinitb), llhoodpij , ys = ysb, pijN = pijsb, sampindex = 1:length(pijsb))

	bhatb <- pijparb$par[2]
	ahatb <- pijparb$par[1]
	phiijhatb <- pijparb$par[3]

	coefsampb <- fixef(glmersampb)
	ranefsampb <- ranef(glmersampb)[[1]][,1]
	names(ranefsampb) <- 1:D

	sigma2ub <- VarCorr(glmersampb)[[1]][1,1]
	checkmessageb <- glmersampb@optinfo$conv$lme4$messages
	check2b <- FALSE
	if(!is.null(checkmessageb)){
	check2b <- checkmessageb == "boundary (singular) fit: see ?isSingular"
	if(check2b){
		itertempb <- 0
		sigma2utempsb <- c()
	repeat{
		itertempb <- itertempb  + 1
		ygentempb <- simulate(glmersampb)
		glmertempb <- glmer(as.matrix(ygentempb )~ xs +  (1 | areafacsamp), family = binomial(link = "logit"), nAGQ = 20)
		sigma2utempsb <-c(sigma2utempsb, VarCorr(glmertempb)[[1]][1,1])
	if(itertempb == 20){break}
	}
	K <- 0.5*sd(sigma2utempsb)
	sigma2ub <- max(c(K, sigma2ub))
	}
	}

	 
	parb <- c( 	coefsampb, sigma2ub, bhatb )
	
	parbs <- rbind(parbs, parb)

	r <- 0
	
	ybarinfrbs <- c()
	yvarinfrbs <-  c()
	yoddsinfrbs <-  c()

	repeat{

		r <- r + 1
		ugenrb <- sapply(sort(unique(areafacsamp)), genudlogist ,parb[3], parb[c(1,2)], xs, ysb, areafacsamp, 100)
		names(ugenrb) <- sort(unique(areafacsamp))
		ugenlongb <- ugenrb[as.character(areafacpop[which(areafacpop %in% areafacsamp)])]
		etapopInfb <- as.vector(cbind(1,xNareasamp)%*%parb[c(1,2)] + ugenlongb - parb[4])
		phatinfb <- exp(etapopInfb)/(1 + exp(etapopInfb))
		ygeninfb <- sapply(phatinfb, function(x){ rbinom(1, prob = x, size = 1)})
		ygeninf[sampindex] <- ysb

		ybarinfrb <- tapply(ygeninfb, areafacpop[which(areafacpop %in% areafacsamp)], mean)
		yvarinfrb <- tapply(ygeninfb, areafacpop[which(areafacpop %in% areafacsamp)], var)
		yoddsinfrb <- tapply(1-ygeninfb, areafacpop[which(areafacpop %in% areafacsamp)], mean)/tapply( ygeninfb, areafacpop[which(areafacpop %in% areafacsamp)], mean)
		
		ybarinfrbs <- rbind(ybarinfrbs, ybarinfrb)
		yvarinfrbs <- rbind(yvarinfrbs, yvarinfrb)
		yoddsinfrbs <- rbind(yoddsinfrbs, yoddsinfrb)

	if(r == R){break}
	}

	ybarinfb <- apply(ybarinfrbs, 2, mean)
	g1infbarb <- apply(ybarinfrbs, 2, var)

	yvarinfb <- apply(yvarinfrbs, 2, mean)
	g1infvarb <- apply(yvarinfrbs, 2, var)

	yoddsinfb <- apply(yoddsinfrbs, 2, mean)
	g1infoddsb <- apply(yoddsinfrbs, 2, var)

      ybarinfbs <- rbind(ybarinfbs, ybarinfb)
	g1infbarbs <- rbind(g1infbarbs, g1infbarb)

	yvarinfbs <- rbind(yvarinfbs, yvarinfb)
	g1infvarbs <- rbind(g1infvarbs, g1infvarb)

	yoddsinfbs <- rbind(yoddsinfbs, yoddsinfb)
	g1infoddsbs <- rbind(g1infoddsbs, g1infoddsb)


	ybarpopbs <-rbind(ybarpopbs,  tapply(yNb, areafacpop, mean))
	yvarpopbs <- rbind(yvarpopbs, tapply(yNb, areafacpop, var))
	yoddspopbs <- rbind(yoddspopbs, tapply(1-yNb, areafacpop, mean)/tapply( yNb, areafacpop, mean)   )

print(paste("boot iter",b))
if(b == B){break}

}


M2hatbar <- apply( ( ybarinfbs  - ybarpopbs)^2, 2, mean)
g1hatbar <- apply(g1infbarbs, 2, mean)

M2hatvar <- apply( (yvarinfbs  - yvarpopbs)^2, 2, mean)
g1hatvar <- apply(g1infvarbs, 2, mean)

M2hatodds <- apply( ( yoddsinfbs  - yoddspopbs)^2, 2, mean)
g1hatodds <- apply(g1infoddsbs, 2, mean)

M2hatbars <- rbind(M2hatbars, M2hatbar)
g1hatbars <- rbind(g1hatbars, g1hatbar)

M2hatvars <- rbind(M2hatvars, M2hatvar)
g1hatvars <- rbind(g1hatvars, g1hatvar)

M2hatoddss <- rbind(M2hatoddss, M2hatodds)
g1hatoddss <- rbind(g1hatoddss, g1hatodds)

vhatpar <- diag(cov(parbs))

vhatpars <- rbind(vhatpars, vhatpar)






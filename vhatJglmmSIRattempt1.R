d <- 0
pards <- c()
repeat{

	d <- d + 1
	ysd <- ys[-which(areafacsamp == d)]
	xsd <- xs[-which(areafacsamp == d)]
	pisd <- piN[sampindex][-which(areafacsamp== d)]
	aread <- areafacsamp[-which(areafacsamp == d)]
	glmersampd <- glmer(cbind(ysd, 1 - ysd) ~ xsd + (1 | aread ), family = binomial(link = "logit"), nAGQ = 20)
	checkmessaged <- glmersampd@optinfo$conv$lme4$messages
	sigma2ud <- VarCorr(glmersampd)[[1]][1,1]
	if(!is.null(checkmessaged)){
	check2d <- checkmessaged  %in% c( "boundary (singular) fit: see ?isSingular", "boundary (singular) fit: see help('isSingular')")
	if(check2d){
		itertempd <- 0
		sigma2utempds <- c()
	repeat{
		itertempd <- itertempd  + 1
		ygentempd <- simulate(glmersampd)
		glmertempd <- glmer(as.matrix(ygentempd )~ xsd +  (1 | aread), family = binomial(link = "logit"), nAGQ = 20)
		sigma2utempds <-c(sigma2utempds, VarCorr(glmertempd)[[1]][1,1])
		if(itertempd == 20){break}
	}
	K <- 0.5*sd(sigma2utempds)
	sigma2ud <- max(c(K, sigma2ud))
	}
	}
	bhat <- lm(log(1/pisd)~xsd + ysd + as.factor(aread) )$coef[3]
	pards <- rbind(pards, c(fixef(glmersampd), transfun(sigma2ud), bhat))
	if(d == D){break}
}

vhatJ <-  cov(pards)*(D-1)^2/D
svdvhatJ <- eigen(vhatJ)
sqrtvhatJ <- svdvhatJ$vectors%*%diag(sqrt(svdvhatJ$values))%*%t(svdvhatJ$vectors)

vhatpars <- rbind(vhatpars, diag(vhatJ))

B <- 100
b <- 0

ybarinfbs <- c()
g1infbarbs <- c()

yvarinfbs <- c()
g1infvarbs <-c()

yoddsinfbs <- c()
g1infoddsbs <-c()


repeat{

	b <- b + 1
      ranvec <- rnorm(ncol(pards), 0, 1)
      ranvec[ranvec > 3] <- 3
	ranvec[ranvec < -3] <- -3
	parb <- c(betahat, transfun(sigma2u), bhat) + as.vector(sqrtvhatJ%*%ranvec)
	parb[3] <- invtransfun(parb[3])
	 
	r <- 0
	
	ybarinfrbs <- c()
	yvarinfrbs <-  c()
	yoddsinfrbs <-  c()

	repeat{

		r <- r + 1
		ugenrb <- sapply(1:D, genudlogist ,parb[3], parb[c(1,2)], xs, ys, areafacsamp, 100)
		names(ugenrb) <- 1:D
		ugenlongb <- ugenrb[as.character(areafacpop)]
		etapopInfb <- as.vector(cbind(1,xN)%*%parb[c(1,2)] + ugenlongb + parb[4])
		phatinfb <- exp(etapopInfb)/(1 + exp(etapopInfb))
		ygeninfb <- sapply(phatinfb, function(x){ rbinom(1, prob = x, size = 1)})
		ybarinfrb <- tapply(ygeninfb, areafacpop, mean)
		yvarinfrb <- tapply(ygeninfb, areafacpop, var)
		yoddsinfrb <- log(tapply(ygeninfb, areafacpop, mean)/tapply(1-ygeninfb, areafacpop, mean))
		
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

print(paste("boot iter",b))

if(b == B){break}

}

M2hatbar <- apply( (t(ybarinfbs) - ybarinf)^2, 1, mean)
g1hatbar <- apply(g1infbarbs, 2, mean)

M2hatvar <- apply( (t(yvarinfbs) - yvarinf)^2, 1, mean)
g1hatvar <- apply(g1infvarbs, 2, mean)

M2hatodds <- apply( (t(yoddsinfbs) - yoddsinf)^2, 1, mean)
g1hatodds <- apply(g1infoddsbs, 2, mean)

M2hatbars <- rbind(M2hatbars, M2hatbar)
g1hatbars <- rbind(g1hatbars, g1hatbar)

M2hatvars <- rbind(M2hatvars, M2hatvar)
g1hatvars <- rbind(g1hatvars, g1hatvar)

M2hatoddss <- rbind(M2hatoddss, M2hatodds)
g1hatoddss <- rbind(g1hatoddss, g1hatodds)














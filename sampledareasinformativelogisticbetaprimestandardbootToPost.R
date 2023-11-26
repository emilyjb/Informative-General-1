
rm(list = ls(all = TRUE))

library("statmod")
library("sampling")
library("rstan")

#library("extraDistr")

source("genudfuns.R")
source("areaparestfuns.R")
source("cdfnonsampu.R")
source("pijparoptfun.R")

D <- 50

Nis <- rep(200, D)

N <- sum(Nis)

nis <- 0.05*Nis
names(nis) <- 1:D

xN <- rnorm(N, mean = 0, sd =1)

areafacpop <- rep(1:D, Nis)

beta0 <- -0.2; beta1 <- 0.7


pbarnoninfs <- c()
pbarinfs <- c()
ybarNis <- c()

iter <- 0 

ybarNis <- c()
yvarNis <- c()
yoddsNis <- c()

ybarnis <- c()
yvarnis <-  c()
yoddsnis <-  c()


ybarinfs <-  c()
ybarnoninfs <-  c()
g1infbars <-  c()

yvarinfs <-  c()
yvarnoninfs <-  c()
g1infvars <- c()

yoddsinfs <- c()
yoddsnoninfs <-  c()
g1infoddss <-  c()


M2hatbars <- c()
g1hatbars <- c()

M2hatvars <- c()
g1hatvars <- c()

M2hatoddss <- c()
g1hatoddss <- c()

areapars <- c()

IDs <- c()

Inftyiter <- c()

unitpars <- c()
vhatpars <- c()

ybpls <- c()
yvpls <- c()
yopls <-  c()

ybplvhats <-  c()
yvplvhats <-  c()
yoplvhats <-  c()

llhoodpijMod <-  function(parinit, ys, pijs, sampindex){
        mupij <- exp(parinit[1] + parinit[2]*ys)/(1 + exp(parinit[1] + parinit[2]*ys))
        shape1 <- mupij*parinit[3] + 1
        shape2 <- (1- mupij)*parinit[3]
        -sum(log( dbeta( pijs,  shape1, shape2)))
}



repeat{

iter <- iter + 1

library("lme4")

 uD <- rnorm(D, mean = 0, sd = 0.5)

names(uD) <- 1:D

uDlong <- uD[as.character(areafacpop)]

etaN <- beta0 + beta1*xN + uDlong

pN <- exp(etaN)/(1 + exp(etaN))

yN <- sapply(pN, function(x){ rbinom(1, prob = x, size = 1)})

gamma2 <- 0.25
pijmean <- exp(-3 + gamma2*yN)/(1 + exp(-3 + gamma2*yN))
phi <- 100
pijN <- rbeta(length(pijmean), shape1 = phi*pijmean, shape2 = phi*(1 - pijmean))
#pinum <- exp( gamma2*yN + 0.5*rnorm(N, mean = 0, sd = 1))
#piden <- tapply(pinum, areafacpop, sum)
#pidenlong <- piden[as.character(areafacpop)]
#nislong <- nis[as.character(areafacpop)]
#piN <- nislong*pinum/pidenlong

IN <- unlist(lapply(as.list(1:D), function(x){ UPsystematic(pijN[areafacpop == x], eps = 1e-18) }))

sampindex <- which(IN == 1)

 
areafacsamp <- areafacpop[sampindex]

ys <- yN[sampindex]
xs <- xN[sampindex]

glmersamp <- glmer(cbind(ys, 1 - ys) ~ xs + (1 | areafacsamp), family = binomial(link = "logit"), nAGQ = 20)
lmwij  <- lm(log(1/pijN[sampindex] -1 )~  ys )
bhat <- lmwij$coef[2]
ahat <- lmwij$coef[1]
mupijinit <- exp(-ahat - bhat*ys)/(1 + exp(-ahat - bhat*ys))
phijinit <- 1/mean((pijN[sampindex] - mupijinit)^2/(mupijinit*(1 - mupijinit)))

pijpar <- optim( c(-ahat, -bhat, phijinit), llhoodpij , ys = ys, pijN = pijN, sampindex = sampindex)

bhat <- pijpar$par[2]
ahat <- pijpar$par[1]
phiijhat <- pijpar$par[3]

coefsamp <- fixef(glmersamp)
ranefsamp <- ranef(glmersamp)[[1]][,1]
names(ranefsamp) <- 1:D

sigma2u <- VarCorr(glmersamp)[[1]][1,1]
checkmessage <- glmersamp@optinfo$conv$lme4$messages
check2 <- FALSE
if(!is.null(checkmessage)){
check2 <- checkmessage %in% c( "boundary (singular) fit: see ?isSingular", "boundary (singular) fit: see help('isSingular')")
if(check2){
itertemp <- 0
sigma2utemps <- c()
repeat{
	itertemp <- itertemp  + 1
	ygentemp <- simulate(glmersamp)
	glmertemp <- glmer(as.matrix(ygentemp )~ xs +  (1 | areafacsamp), family = binomial(link = "logit"), nAGQ = 20)
	sigma2utemps <-c(sigma2utemps, VarCorr(glmertemp)[[1]][1,1])
	if(itertemp == 20){break}
}
K <- 0.5*sd(sigma2utemps)
sigma2u <- max(c(K, sigma2u))
}
}

XN <- cbind(1, xN)
 
ybarNis <- rbind(ybarNis, tapply(yN, areafacpop, mean))
yvarNis <- rbind(yvarNis, tapply(yN, areafacpop, var))
#yoddsNis <- rbind(yoddsNis, log( tapply(yN, areafacpop, mean)/tapply(1-yN, areafacpop, mean) ) )
yoddsNis <- rbind(yoddsNis,   tapply(1-yN, areafacpop, mean)/tapply( yN, areafacpop, mean) ) 

ybarnis <- rbind(ybarnis, tapply(ys, areafacsamp, mean))
yvarnis <- rbind(yvarnis, tapply(ys, areafacsamp, var))
#yoddsnis <- rbind(yoddsnis, log( tapply(ys, areafacsamp, mean)/tapply(1-ys, areafacsamp, mean) ) )
yoddsnis <- rbind(yoddsnis,  tapply(1-ys, areafacsamp, mean)/tapply( ys, areafacsamp, mean) ) 

unitpars <- rbind(unitpars, c(coefsamp, sigma2u, bhat))

#source("informativelogisticpredranSIR.R")
#source("informativelogisticpredranSIRSAreasBetaPrime.R")
source("informativelogisticpredranSIRSampledAreasBetaPrime.R")

#source("vhatJglmmSIRattempt1.R")
#source("vhatJglmmSIRattempt1NSAreas.R")

source("sampledareasbootstrapstandardlogisticbetaprime.R")
#source("logisticbootstrapbetaprime.R")
#source("bootstrapstandardlogisticbetaprime.R")

#IDs <- rbind(IDs, ID)

piN <- pijN
source("runpl.R")

print(paste(iter))

if(iter%%10 == 0){ save.image("betaprimestandardboot11-5-23.Rdata")}

if(iter == 500){break}


}















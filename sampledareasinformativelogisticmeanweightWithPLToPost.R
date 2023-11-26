
rm(list = ls(all = TRUE))

library("statmod")
library("sampling")
library("rstan")

source("genudfuns.R")
source("areaparestfuns.R")
source("cdfnonsampu.R")

### Comment out : D <- 50
D <- 50

Nis <- rep(200, D)

N <- sum(Nis)

nis <- 0.05*Nis
names(nis) <- 1:D

xN <- rnorm(N, mean = 0, sd =1)

areafacpop <- rep(1:D, Nis)

beta0 <- -0.2 ; beta1 <- 0.7


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
pinum <- exp( gamma2*yN + 0.5*rnorm(N, mean = 0, sd = 1))
piden <- tapply(pinum, areafacpop, sum)
pidenlong <- piden[as.character(areafacpop)]
nislong <- nis[as.character(areafacpop)]
piN <- nislong*pinum/pidenlong

IN <- unlist(lapply(as.list(1:D), function(x){ UPsystematic(piN[areafacpop == x], eps = 1e-18) }))

sampindex <- which(IN == 1)

areafacsamp <- areafacpop[sampindex]

ys <- yN[sampindex]
xs <- xN[sampindex]

glmersamp <- glmer(cbind(ys, 1 - ys) ~ xs + (1 | areafacsamp), family = binomial(link = "logit"), nAGQ = 20)
bhat <- lm(log(1/piN[sampindex])~xs + ys + as.factor(areafacsamp) )$coef[3]

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

ws <- 1/piN[sampindex]
ybarnis <- rbind(ybarnis, tapply(ys*ws, areafacsamp, sum)/tapply(ws, areafacsamp, sum))
yvarnis <- rbind(yvarnis, tapply(ys^2*ws, areafacsamp, sum)/tapply(ws, areafacsamp, sum) - (tapply(ys*ws, areafacsamp, sum)/tapply(ws, areafacsamp, sum))^2    )
#yoddsnis <- rbind(yoddsnis, log( tapply(ys, areafacsamp, mean)/tapply(1-ys, areafacsamp, mean) ) )
yoddsnis <- rbind(yoddsnis,  tapply( (1-ys)*ws, areafacsamp, mean)/tapply( ys*ws, areafacsamp, mean) ) 

unitpars <- rbind(unitpars, c(coefsamp, sigma2u, bhat))

source("informativelogisticpredranSIR.R")
#source("informativelogisticpredranSIRNSAreas.R")

source("vhatJglmmSIRattempt1.R")
#source("vhatJglmmSIRattempt1NSAreas.R")
#source("vhatJglmmSIRattempt2NSAreasME.R")

#IDs <- rbind(IDs, ID)



source("runpl.R")


if(iter%%100 == 0){save.image("meanweight11-12-23.Rdata")}

print(paste(iter))

if(iter == 500){break}


}
















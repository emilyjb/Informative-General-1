
ws <- 1/piN[sampindex]
wtildes <- length(ws)*ws/sum(ws)

datalist <- list(D = D, N = length(ws), P = 2, Y = ys, wgt = wtildes, cty = areafacsamp, X = cbind(1, xs))

pseudofit <- stan(file = "pseudologistparker.stan", data = datalist, chains = 3, iter = 1100, warmup = 100)

parameters <- extract(pseudofit)

iterpost <- 0
ybps <- c()
yvps <- c()
yops <-  c()

repeat{

iterpost <- iterpost + 1

betahatb <- parameters$beta[iterpost,]
uhatb <- parameters$mu[iterpost,]

names(uhatb) <- 1:D

ran <- uhatb[as.character(areafacpop)]
fix <- (XN%*%betahatb)[,1]

linkpop <- fix + ran
mupop <- exp(linkpop)/(1 + exp(linkpop))

ypop <- sapply(mupop, function(x){ rbinom(1, prob = x, size = 1)})
ypop[sampindex] <- ys

ybp <- tapply(ypop, areafacpop, mean)
yvp <- tapply(ypop, areafacpop, var)
yop <- tapply(1- ypop, areafacpop, mean)/tapply(ypop, areafacpop, mean)

ybps <-  rbind(ybps, ybp)
yvps <- rbind(yvps, yvp)
yops <- rbind(yops, yop)

if(iterpost == nrow(parameters$beta)){break}

}


ybpl <- apply(ybps, 2, mean)
yvpl <- apply(yvps, 2, mean)
yopl <- apply(yops, 2, mean)

ybplvhat <- apply(ybps, 2, var)
yvplvhat <- apply(yvps, 2, var)
yoplvhat <- apply(yops, 2, var)


ybpls <- rbind(ybpls, ybpl)
yvpls <- rbind(yvpls, yvpl)
yopls <- rbind(yopls, yopl)

ybplvhats <- rbind(ybplvhats, ybplvhat)
yvplvhats <- rbind(yvplvhats, yvplvhat)
yoplvhats <- rbind(yoplvhats, yoplvhat)


















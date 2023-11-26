

maxcnt <- 500

r1 <- apply(unitpars, 2, mean)
r2 <- apply(unitpars, 2, var)  
r3 <- apply(vhatpars, 2, mean)
 
round(rbind(r1, r2, r3), digits = 3)

rbinfbar <- mean(abs(apply( (ybarinfs[1:maxcnt,] - ybarNis[1:maxcnt,]) , 2, mean)/apply( ( ybarNis[1:maxcnt,]) , 2, mean) ))
rbnoninfbar <-   mean(abs(apply( (ybarnoninfs[1:maxcnt,] - ybarNis[1:maxcnt,]) , 2, mean)/apply( ( ybarNis[1:maxcnt,]) , 2, mean) ))

msemultbar <- g1infbars[1:maxcnt,]^2/g1hatbars[1:maxcnt,] + M2hatbars[1:maxcnt,]
mseaddbar <-  2*g1infbars[1:maxcnt,] - g1hatbars[1:maxcnt,] + M2hatbars[1:maxcnt,]
msecompbar <- ifelse(2*g1infbars[1:maxcnt,] - g1hatbars[1:maxcnt,] > 0, mseaddbar[1:maxcnt,], msemultbar[1:maxcnt,])
msenbar <- g1infbars[1:maxcnt,] + M2hatbars[1:maxcnt,]
msenoninfbar <- mean(apply( (ybarnoninfs[1:maxcnt,] - ybarNis[1:maxcnt,])^2, 2, mean))
mseinfbar <- mean(apply( (ybarinfs[1:maxcnt,] - ybarNis[1:maxcnt,])^2, 2, mean))
msehatbar <- mean(M2hatbars[1:maxcnt,])

rbnoninfbar
rbinfbar
msenoninfbar
mseinfbar
msehatbar

rbinfvar <- mean(abs(apply( (yvarinfs[1:maxcnt,] - yvarNis[1:maxcnt,]) , 2, mean)/apply( ( yvarNis[1:maxcnt,]) , 2, mean) ))
rbnoninfvar <-   mean(abs(apply( (yvarnoninfs[1:maxcnt,] - yvarNis[1:maxcnt,]) , 2, mean)/apply( ( yvarNis[1:maxcnt,]) , 2, mean) ))
msemultvar <- g1infvars[1:maxcnt,]^2/g1hatvars[1:maxcnt,] + M2hatvars[1:maxcnt,]
mseaddvar <-  2*g1infvars[1:maxcnt,] - g1hatvars[1:maxcnt,] + M2hatvars[1:maxcnt,]
msecompvar <- ifelse(2*g1infvars[1:maxcnt,] - g1hatvars[1:maxcnt,] > 0, mseaddvar[1:maxcnt,], msemultvar[1:maxcnt,])
msenvar <- g1infvars[1:maxcnt,] + M2hatvars[1:maxcnt,]
msenoninfvar <- mean(apply( (yvarnoninfs[1:maxcnt,] - yvarNis[1:maxcnt,])^2, 2, mean))
mseinfvar <- mean(apply( (yvarinfs[1:maxcnt,] - yvarNis[1:maxcnt,])^2, 2, mean))
msehatvar <- mean(M2hatvars[1:maxcnt,])

rbnoninfvar
rbinfvar
msenoninfvar
mseinfvar
msehatvar



rbinfodds <- mean(abs(apply( (yoddsinfs[1:maxcnt,] - yoddsNis[1:maxcnt,]) , 2, mean)/apply( ( yoddsNis[1:maxcnt,]) , 2, mean) ))
rbnoninfodds <-   mean(abs(apply( (yoddsnoninfs[1:maxcnt,] - yoddsNis[1:maxcnt,]) , 2, mean)/apply( ( yoddsNis[1:maxcnt,]) , 2, mean) ))
msemultodds <- g1infoddss[1:maxcnt,]^2/g1hatoddss[1:maxcnt,] + M2hatoddss[1:maxcnt,]
mseaddodds <-  2*g1infoddss[1:maxcnt,] - g1hatoddss[1:maxcnt,] + M2hatoddss[1:maxcnt,]
msecompodds <- ifelse(2*g1infoddss[1:maxcnt,] - g1hatoddss[1:maxcnt,] > 0, mseaddodds[1:maxcnt,], msemultodds[1:maxcnt,])
msenodds <- g1infoddss[1:maxcnt,] + M2hatoddss[1:maxcnt,]
msenoninfodds <- mean(apply( (yoddsnoninfs[1:maxcnt,] - yoddsNis[1:maxcnt,])^2, 2, mean))
mseinfodds <- mean(apply( (yoddsinfs[1:maxcnt,] - yoddsNis[1:maxcnt,])^2, 2, mean))
msehatodds <- mean( M2hatoddss[1:maxcnt,])

rbnoninfodds
rbinfodds
msenoninfodds
mseinfodds
msehatodds


c1 <- matrix(c(rbnoninfbar, rbinfbar, rbnoninfvar, rbinfvar, rbnoninfodds, rbinfodds), ncol = 2, byrow= TRUE)
c2 <- matrix(c(msenoninfbar,mseinfbar,msehatbar,msenoninfvar,mseinfvar,msehatvar,msenoninfodds,mseinfodds,msehatodds), ncol = 3, byrow = TRUE)

library("xtable")
xtable(round(100*cbind(c1, c2), digits = 3))


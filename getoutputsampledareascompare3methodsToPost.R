

maxcnt <- 500

# Output for the mean 

b1m <- apply( (ybarinfs[1:maxcnt,] - ybarNis[1:maxcnt,]) , 2, mean)/apply(ybarNis[1:maxcnt,], 2, mean)
b2m <- apply( (ybarnoninfs[1:maxcnt,] - ybarNis[1:maxcnt,]) , 2, mean)/apply(ybarNis[1:maxcnt,], 2, mean)
b3m <- apply( (ybpls[1:maxcnt,] - ybarNis[1:maxcnt,]), 2, mean)/apply(ybarNis[1:maxcnt,], 2, mean) 

boxplot(b1m, b2m, b3m)
abline(h = 0)

ab1m <- mean(abs(apply( (ybarinfs[1:maxcnt,] - ybarNis[1:maxcnt,]) , 2, mean)/apply( ( ybarNis[1:maxcnt,]) , 2, mean) ))
ab2m <- mean(abs(apply( (ybarnoninfs[1:maxcnt,] - ybarNis[1:maxcnt,]) , 2, mean)/apply( ( ybarNis[1:maxcnt,]) , 2, mean) ))
ab3m <- mean(abs(apply( (ybpls[1:maxcnt,] - ybarNis[1:maxcnt,]) , 2, mean)/apply( ( ybarNis[1:maxcnt,]) , 2, mean) ))

mse1m <- mean(apply( (ybarinfs[1:maxcnt,] - ybarNis[1:maxcnt,])^2, 2, mean))
mse2m <- mean(apply( (ybarnoninfs[1:maxcnt,] - ybarNis[1:maxcnt,])^2, 2, mean))
mse3m <- mean(apply( (ybpls[1:maxcnt,] - ybarNis[1:maxcnt,])^2, 2, mean))

msemultbar <- g1infbars[1:maxcnt,]^2/g1hatbars[1:maxcnt,] + M2hatbars[1:maxcnt,]
mseaddbar <-  2*g1infbars[1:maxcnt,] - g1hatbars[1:maxcnt,] + M2hatbars[1:maxcnt,]
msecompbar <- ifelse(mseaddbar > 0, mseaddbar, msemultbar)

rb1bar <- 100*(apply(msecompbar, 2, mean)/apply( (ybarinfs[1:maxcnt,] - ybarNis[1:maxcnt,])^2, 2, mean)- 1)
rb3bar <- 100*(apply(ybplvhats, 2, mean)/apply( (ybpls[1:maxcnt,] - ybarNis[1:maxcnt,])^2, 2, mean) - 1)

100*(mean(msecompbar)/mse1m - 1)
100*(mean(ybplvhats)/mse3m - 1)

 

# Output for the variance

b1v <- apply( (yvarinfs[1:maxcnt,] - yvarNis[1:maxcnt,]) , 2, mean)/apply(yvarNis[1:maxcnt,], 2, mean)
b2v <- apply( (yvarnoninfs[1:maxcnt,] - yvarNis[1:maxcnt,]) , 2, mean)/apply(yvarNis[1:maxcnt,], 2, mean)
b3v <- apply( (yvpls[1:maxcnt,] - yvarNis[1:maxcnt,]), 2, mean)/apply(yvarNis[1:maxcnt,], 2, mean) 

boxplot(b1v, b2v, b3v)
abline(h = 0)

ab1v <- mean(abs(apply( (yvarinfs[1:maxcnt,] - yvarNis[1:maxcnt,]) , 2, mean)/apply( ( yvarNis[1:maxcnt,]) , 2, mean) ))
ab2v <- mean(abs(apply( (yvarnoninfs[1:maxcnt,] - yvarNis[1:maxcnt,]) , 2, mean)/apply( ( yvarNis[1:maxcnt,]) , 2, mean) ))
ab3v <- mean(abs(apply( (yvpls[1:maxcnt,] - yvarNis[1:maxcnt,]) , 2, mean)/apply( ( yvarNis[1:maxcnt,]) , 2, mean) ))

mse1v <- mean(apply( (yvarinfs[1:maxcnt,] - yvarNis[1:maxcnt,])^2, 2, mean))
mse2v <- mean(apply( (yvarnoninfs[1:maxcnt,] - yvarNis[1:maxcnt,])^2, 2, mean))
mse3v <- mean(apply( (yvpls[1:maxcnt,] - yvarNis[1:maxcnt,])^2, 2, mean))

msemultvar <- g1infvars[1:maxcnt,]^2/g1hatvars[1:maxcnt,] + M2hatvars[1:maxcnt,]
mseaddvar <-  2*g1infvars[1:maxcnt,] - g1hatvars[1:maxcnt,] + M2hatvars[1:maxcnt,]
msecompvar <- ifelse(mseaddvar > 0, mseaddvar, msemultvar)

rb1var <- 100*(apply(msecompvar, 2, mean)/apply( (yvarinfs[1:maxcnt,] - yvarNis[1:maxcnt,])^2, 2, mean)- 1)
rb3var <- 100*(apply(yvplvhats, 2, mean)/apply( (yvpls[1:maxcnt,] - yvarNis[1:maxcnt,])^2, 2, mean) - 1)

100*(mean(msecompvar)/mse1v - 1)
100*(mean(yvplvhats)/mse3v - 1)
 
# Output for the odds ratio: 

b1o <- apply( (yoddsinfs[1:maxcnt,] - yoddsNis[1:maxcnt,]) , 2, mean)/apply(yoddsNis[1:maxcnt,], 2, mean)
b2o <- apply( (yoddsnoninfs[1:maxcnt,] - yoddsNis[1:maxcnt,]) , 2, mean)/apply(yoddsNis[1:maxcnt,], 2, mean)
b3o <- apply( (yopls[1:maxcnt,] - yoddsNis[1:maxcnt,]), 2, mean)/apply(yoddsNis[1:maxcnt,], 2, mean) 

boxplot(b1o, b2o, b3o)
abline(h = 0)

ab1o <- mean(abs(apply( (yoddsinfs[1:maxcnt,] - yoddsNis[1:maxcnt,]) , 2, mean)/apply( ( yoddsNis[1:maxcnt,]) , 2, mean) ))
ab2o <- mean(abs(apply( (yoddsnoninfs[1:maxcnt,] - yoddsNis[1:maxcnt,]) , 2, mean)/apply( ( yoddsNis[1:maxcnt,]) , 2, mean) ))
ab3o <- mean(abs(apply( (yopls[1:maxcnt,] - yoddsNis[1:maxcnt,]) , 2, mean)/apply( ( yoddsNis[1:maxcnt,]) , 2, mean) ))

mse1o <- mean(apply( (yoddsinfs[1:maxcnt,] - yoddsNis[1:maxcnt,])^2, 2, mean))
mse2o <- mean(apply( (yoddsnoninfs[1:maxcnt,] - yoddsNis[1:maxcnt,])^2, 2, mean))
mse3o <- mean(apply( (yopls[1:maxcnt,] - yoddsNis[1:maxcnt,])^2, 2, mean))

msemultodds <- g1infoddss[1:maxcnt,]^2/g1hatoddss[1:maxcnt,] #+ M2hatoddss[1:maxcnt,]
mseaddodds <-  2*g1infoddss[1:maxcnt,] - g1hatoddss[1:maxcnt,]# + M2hatoddss[1:maxcnt,]
msecompodds <- ifelse(mseaddodds > 0, mseaddodds, msemultodds)

rb1odds <- 100*(apply(g1infoddss[1:maxcnt,], 2, mean)/apply( (yoddsinfs[1:maxcnt,] - yoddsNis[1:maxcnt,])^2, 2, mean)- 1)
rb3odds <- 100*(apply(yoplvhats, 2, mean)/apply( (yopls[1:maxcnt,] - yoddsNis[1:maxcnt,])^2, 2, mean) - 1)

100*(mean(g1infoddss[1:maxcnt,])/mse1o - 1)
100*(mean(yoplvhats)/mse3o - 1)

boxplot(cbind(rb1odds, rb3odds))


dfbox <- data.frame(b1m, b2m,b3m,  b1v, b2v,b3v,  b1o, b2o, b3o)
colnames(dfbox) <- c("Mean Inf", "Mean Non-Inf", "Mean PL", "Var Inf", "Var Non-Inf", "Var PL", "Odd Inf", "Odds Non-Inf", "Odds PL")

boxplot(dfbox)


rm <- c(ab1m, ab2m, ab3m, mse1m, mse2m, mse3m, mean(msecompbar), mean(ybplvhats))
rv <- c(ab1v, ab2v, ab3v, mse1v, mse2v, mse3v, mean(msecompvar), mean(yvplvhats))
ro <- c(ab1o, ab2o, ab3o, mse1o, mse2o, mse3o, mean(g1infoddss), mean(yoplvhats))

library("xtable")
xtable(round(rbind(rm, rv, ro)*100, digits = 2), digits = 2)










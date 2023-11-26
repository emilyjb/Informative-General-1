


llhoodpij <- function(parinit, ys, pijN, sampindex){
	mupij <- exp(parinit[1] + parinit[2]*ys)/(1 + exp(parinit[1] + parinit[2]*ys))
	shape1 <- mupij*parinit[3] + 1
	shape2 <- (1- mupij)*parinit[3]
	-sum(log( dbeta( pijN[sampindex],  shape1, shape2)))
}




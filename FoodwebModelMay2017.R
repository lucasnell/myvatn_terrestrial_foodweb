##Food-web model for Myvatn (May 2017)

# This is a nutrient-balance model which was much more dynamically stable than a hybrid nutrient/population model that I tried first. There is mass balance, with nutrients "leaking" from the system through leaching and entering the system through external inputs and midges that enter the detritus pool.

# Note that this is an example in which a nutrient balance model is inherently more dyanmically stable than a population-level model. This topic (relative stability of mass-balance vs. population models) is something Joe and I talked about looking into a long time ago.

#install.packages("rootSolve")
library(rootSolve)

#########################################################
# functions
#########################################################
# all 6 variates
tosolvevector <- function(par, pmatrix){
	pmatrix[is.na(pmatrix)] <- par
	a <- pmatrix[,1:6]
	x <- pmatrix[,7]
	b <- pmatrix[,8]
	N <- b[3]
	M <- b[4]
	
	s1 <- -a[1,1]*x[1] + a[1,2]*x[2]*x[3] - a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) + N
	s2 <- -a[1,2]*x[2]*x[3] - a[3,2]*x[2]*x[3] + a[3,3]*x[3] + a[4,4]*x[4] + a[5,5]*x[5] + a[6,6]*x[6] + M
	s3 <- a[3,2]*x[2] - a[3,6]*x[6]*(1 + b[2] * x[6])^(-1) - a[3,3]
	s4 <- a[1,4]*x[1]*(1 + b[1]*x[4])^(-1) - a[4,5]*x[5] - a[4,4]
	s5 <- a[4,5]*x[4] - a[5,6]*x[6]*(1 + b[2] * x[6])^(-1) - a[5,5]
	s6 <- a[3,6]*x[3]*(1 + b[2] * x[6])^(-1) + a[5,6]*x[5]*(1 + b[2] * x[6])^(-1) - a[6,6]

	return(c(s1,s2,s3,s4,s5,s6)) 
}

# removing herbivores and predators
tosolvevector4 <- function(par, pmatrix){
	pmatrix[is.na(pmatrix)] <- par
	a <- pmatrix[,1:4]
	x <- pmatrix[,5]
	b <- pmatrix[,6]
	N <- b[3]
	M <- b[4]
	
	s1 <- -a[1,1]*x[1] + a[1,2]*x[2]*x[3] - a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) + N
	s2 <- -a[1,2]*x[2]*x[3] - a[3,2]*x[2]*x[3] + a[3,3]*x[3] + a[4,4]*x[4] + M
	s3 <- a[3,2]*x[2] - a[3,3]
	s4 <- a[1,4]*x[1]*(1 + b[1]*x[4])^(-1) - a[4,4]

	return(c(s1,s2,s3,s4)) 
}
#########################################################
#########################################################
# all 6 variables
#########################################################
#########################################################

# the code below solves for parameters given x

a <- matrix(0, nrow=6, ncol=6)
x <- matrix(0, nrow=6, ncol=1)
b <- matrix(0, nrow=6, ncol=1)

# 1-nutrients
# 2-detritus
# 3-detritivores
# 4-plants
# 5-herbivores
# 6-predators

# each element a[i,j] gives a transition rate for i from j. The a[i,i] elements are mortality or loss. The equations are set up with the right sign, so all coefficients are positive.

x[1] <- 343000
x[2] <- 114000
x[3] <- 81
x[4] <- 4300
x[5] <- 24
x[6] <- 13

a[1,1] <- 100/x[1]
a[3,3] <- 2/x[3]
a[4,4] <- .1/x[4]
a[5,5] <- 2/x[5]
a[6,6] <- 2/x[6]

a[1,2] <- NA
a[3,2] <- NA
a[1,4] <- NA
a[4,5] <- NA
a[3,6] <- NA
a[5,6] <- NA

b[1] <- 100
b[2] <- .01
#N
b[3] <- 100
#M
b[4] <- 0

start=c(1/(x[2]*x[3]), a[3,3]/x[2], 1/x[1], 1/x[4], 1/x[3], 1/x[5])
pmatrix <- cbind(a,x,b)
tosolvevector(start, pmatrix)
z <- multiroot(f=tosolvevector, start=start, parms = pmatrix, maxiter = 100, rtol = 1e-6, atol = 1e-8, ctol = 1e-8, useFortran = TRUE, positive = T,jacfunc = NULL, jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1)

tosolvevector(start, pmatrix)
# This contains all of the parameters 
pfitted <- pmatrix
pfitted[is.na(pmatrix)] <- z$root

#########################################################
# simulate model at equilibrium with M=0

tscale <- 1
a <- tscale*pfitted[,1:6]
x <- pfitted[,7]
b <- pfitted[,8]
N <- tscale*b[3]
M <- b[4]

nt <- 500
X <- matrix(0, nrow=nt, ncol=6)

for(t in 1:nt){
	x1 <- x[1] - a[1,1]*x[1] + a[1,2]*x[2]*x[3] - a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) + N
	x2 <- x[2] - a[1,2]*x[2]*x[3] - a[3,2]*x[2]*x[3] + a[3,3]*x[3] + a[4,4]*x[4] + a[5,5]*x[5] + a[6,6]*x[6] + M
	x3 <- x[3] + a[3,2]*x[2]*x[3] - a[3,6]*x[3]*x[6]*(1 + b[2] * x[6])^(-1) - a[3,3]*x[3]
	x4 <- x[4] + a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) - a[4,5]*x[4]*x[5] - a[4,4]*x[4]
	x5 <- x[5] + a[4,5]*x[4]*x[5] - a[5,6]*x[5]*x[6]*(1 + b[2] * x[6])^(-1) - a[5,5]*x[5]
	x6 <- x[6] + a[3,6]*x[3]*x[6]*(1 + b[2] * x[6])^(-1) + a[5,6]*x[5]*x[6]*(1 + b[2] * x[6])^(-1) - a[6,6]*x[6]
	
	x <- c(x1, x2, x3, x4, x5, x6)
	X[t,] <- x
}

par(mfrow=c(2,2))
matplot(X, typ="l", lty = 1, xaxt="n", xlab="Time", log='y', ylab="log Abundance", col=1:6)

# This is my rendition of a foodweb picture, with circles proportional to abundance
plot(c(0,-1,-1,1,1,0),c(0,1,2,1,2,3), cex=log10(X[nt,]), axes=F, xlim=c(-2,2), ylim=c(-1,4), xlab="", ylab="", main="log Abundance", col=1:6)

#########################################################
# simulate model for experiment
# This starts at equilibrium and adds midges as a press
startX <- X[nt,]
x <- startX
M <- 0
maxM <- 10

nt <- 20000
X <- matrix(0, nrow=nt, ncol=6)
for(t in 1:nt){
	
	if(t > nt/10) M <- maxM
	
	x1 <- x[1] - a[1,1]*x[1] + a[1,2]*x[2]*x[3] - a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) + N
	x2 <- x[2] - a[1,2]*x[2]*x[3] - a[3,2]*x[2]*x[3] + a[3,3]*x[3] + a[4,4]*x[4] + a[5,5]*x[5] + a[6,6]*x[6] + M
	x3 <- x[3] + a[3,2]*x[2]*x[3] - a[3,6]*x[3]*x[6]*(1 + b[2] * x[6])^(-1) - a[3,3]*x[3]
	x4 <- x[4] + a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) - a[4,5]*x[4]*x[5] - a[4,4]*x[4]
	x5 <- x[5] + a[4,5]*x[4]*x[5] - a[5,6]*x[5]*x[6]*(1 + b[2] * x[6])^(-1) - a[5,5]*x[5]
	x6 <- x[6] + a[3,6]*x[3]*x[6]*(1 + b[2] * x[6])^(-1) + a[5,6]*x[5]*x[6]*(1 + b[2] * x[6])^(-1) - a[6,6]*x[6]
	
	x <- c(x1, x2, x3, x4, x5, x6)
	X[t,] <- x
}

matplot(X, typ="l", lty = 1, xaxt="n", xlab="Time", log='y', ylab="log Abundance", col=1:6)
# This plots the foodweb in terms of log change in abundance
plot(c(0,-1,-1,1,1,0),c(0,1,2,1,2,3), cex=5*log(abs(X[nt,]/startX)), axes=F, xlim=c(-2,2), ylim=c(-1,4), xlab="", ylab="", main="log Change in abundance", col=1:6)




#########################################################
#########################################################
# model without predators and herbivores
#########################################################
#########################################################

# the code below solves for parameters given x

a <- matrix(0, nrow=4, ncol=4)
x <- matrix(0, nrow=4, ncol=1)
b <- matrix(0, nrow=4, ncol=1)

# 1-nutrients
# 2-detritus
# 3-detritivores
# 4-plants

# each element a[i,j] gives a transition rate for i from j. The a[i,i] elements are mortality or loss. The equations are set up with the right sign, so all coefficients are positive.

x[1] <- 343000
x[2] <- 114000
x[3] <- 81
x[4] <- 4300

a[1,1] <- 100/x[1]
a[3,3] <- 2/x[3]

a[1,2] <- NA
a[3,2] <- NA
a[1,4] <- NA
a[4,4] <- NA

b[1] <- 100
#N
b[3] <- 100
#M
b[4] <- 0

start=c(1/(x[2]*x[3]), a[3,3]/x[2], 1/x[1], 1/x[4])

pmatrix <- cbind(a,x,b)
z <- multiroot(f=tosolvevector4, start=start, parms = pmatrix, maxiter = 100, rtol = 1e-6, atol = 1e-8, ctol = 1e-8, useFortran = TRUE, positive = T,jacfunc = NULL, jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1)
z

# This contains all of the parameters 
pfitted <- pmatrix
pfitted[is.na(pmatrix)] <- z$root

#########################################################
# simulate model at equilibrium with M=0

tscale <- 1
a <- tscale*pfitted[,1:4]
x <- pfitted[,5]
b <- pfitted[,6]
N <- tscale*b[3]
M <- b[4]

nt <- 500
X <- matrix(0, nrow=nt, ncol=4)

for(t in 1:nt){
	x1 <- x[1] - a[1,1]*x[1] + a[1,2]*x[2]*x[3] - a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) + N
	x2 <- x[2] - a[1,2]*x[2]*x[3] - a[3,2]*x[2]*x[3] + a[3,3]*x[3] + a[4,4]*x[4] + M
	x3 <- x[3] + a[3,2]*x[2]*x[3] - a[3,3]*x[3]
	x4 <- x[4] + a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) - a[4,4]*x[4]
	
	x <- c(x1, x2, x3, x4)
	X[t,] <- x
}

par(mfrow=c(2,2))
matplot(X, typ="l", lty = 1, xaxt="n", xlab="Time", log='y', ylab="log Abundance", col=1:4)

# This is my rendition of a foodweb picture, with circles proportional to abundance
plot(c(0,-1,-1,1),c(0,1,2,1), cex=log10(X[nt,]), axes=F, xlim=c(-2,2), ylim=c(-1,4), xlab="", ylab="", main="log Abundance", col=1:4)

#########################################################
# simulate model for experiment
# This starts at equilibrium and adds midges as a press
startX <- X[nt,]
x <- startX
M <- 0
maxM <- 10

nt <- 20000
X <- matrix(0, nrow=nt, ncol=4)
for(t in 1:nt){
	
	if(t > nt/10) M <- maxM
	
	x1 <- x[1] - a[1,1]*x[1] + a[1,2]*x[2]*x[3] - a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) + N
	x2 <- x[2] - a[1,2]*x[2]*x[3] - a[3,2]*x[2]*x[3] + a[3,3]*x[3] + a[4,4]*x[4] + M
	x3 <- x[3] + a[3,2]*x[2]*x[3] - a[3,3]*x[3]
	x4 <- x[4] + a[1,4]*x[1]*x[4]*(1 + b[1]*x[4])^(-1) - a[4,4]*x[4]
	
	x <- c(x1, x2, x3, x4)
	X[t,] <- x
}

matplot(X, typ="l", lty = 1, xaxt="n", xlab="Time", log='y', ylab="log Abundance", col=1:4)
# This plots the foodweb in terms of log change in abundance
plot(c(0,-1,-1,1),c(0,1,2,1), cex=2*log(abs(X[nt,]/startX)), axes=F, xlim=c(-2,2), ylim=c(-1,4), xlab="", ylab="", main="log Change in abundance", col=1:4)



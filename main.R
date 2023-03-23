#========================================================== 
# Estimate VECM and compute Information Shares based on ICA
#==========================================================

# Clear environment and console
rm(list = ls(all = TRUE))        # clear environment 
cat("\014")                      # clear console

setwd("C:/Users/sebas/Documents/Università&Ricerca/Ricerca/Ricerca personale/Price Discovery/IC-price discovery/Rcodes")
source("set.of.procedures.R")
source("useful functions.R")

# import data (depending on the price discovery analysis you want to make)
data <- read.csv('./data processed/df.qt.1sec.IBM.csv', sep = ',', header = TRUE)
Y <- as.matrix(log(data[,-c(1)]))  # get log-prices (remove from data the columns not needed)

# VECM estimation
T <- dim(Y)[1]  # time series dimension              
n <- dim(Y)[2]  # cross section dimension
nb.lags = 10    # number of lags in the model
B <- cbind(matrix(c(1,1,1), ncol = 1,  nrow = dim(Y)[2]-1), - diag(dim(Y)[2]-1)) # cointegrating matrix is known 

# lag variables
X <- NULL
for(i in 1:nb.lags){
  lagged.Y <- rbind(
    matrix(NaN,i,n),
    Y[1:(T-i),])
  X <- cbind(X,lagged.Y)
}

# take first differences ('small' stands for 'differentiated variable', 'capital' stands for 'variable in levels' )
y <- diff(Y, lag = 1) 
x <- diff(X, lag = 1)
BY <- t(B%*%t(Y))[-nrow(Y),] # time series of the cointegrating relationships among price series

Phi <- matrix(0,n,n*nb.lags)
a <- matrix(0,n,n-1)
Eta <- NULL # Eta is the matrix of OLS residuals

# estimate equations and get residuals:
for(i in 1:n){
  eq <- lm(y[,i] ~ BY + x + 0)
  Phi[i,] <- eq$coef[n:(dim(Phi)[2]+(n-1))] 
  a[i,] <- eq$coef[1:(n-1)]
  Eta <- cbind(Eta,eq$residuals)
}

# Get ingredients to compute Hasbrouck's IS (Psi will be needed also for the IC-IS measure)
library(mcompanion)

phi_k <- split_mat(Phi,n,n)
phi_sum <- matrix(0,n,n)

for (i in 1:nb.lags) {
  phi_sum <- phi_sum + phi_k[1:n, 1:n, i]
}

a_orth <- null_complement(a)
B_orth <- null_complement(t(B))
I <- diag(n)
H <- solve(t(a_orth)%*%(I-phi_sum)%*%B_orth)
Psi1 <- B_orth%*%H%*%t(a_orth)
Sigma_eta <- cov(Eta) #covariance matrix of OLS residuals



# -------------------------------------------------------------------------
# Estimate C for the IC-IS measure, PML approach of Gourieroux et al. 2017
# -------------------------------------------------------------------------

# decide the distribution to be used in the pseudo maximum likelihood estimation and set useful parameters
distri <- list(
  type=c("hyper.sec","hyper.sec","hyper.sec", "hyper.sec"),
  df=c(3,4,3,4),
  p=c(.5,.5,.5,.5),
  mu=c(.1,.1,.1,.1),
  sigma=c(.5,.7,.9,1.1)
)

step <- .0001
xxxx <- seq(-6,6,by=step)
yyyy <- exp(g(matrix(xxxx,length(xxxx),4),distri,indic.deriv=0)$log.L)
par(mfrow=c(1,1))
plot(xxxx,yyyy[,1],type="l")
lines(xxxx,yyyy[,2],col="red")
lines(xxxx,yyyy[,3],col="blue")
lines(xxxx,yyyy[,4],col="gray")

XXXX <- matrix(xxxx,length(xxxx),4)
XY <- XXXX * yyyy
m1 <- apply(XY,2,function(x){sum(x)*step})
m2 <- apply(XXXX^2 * yyyy,2,function(x){sum(x)*step})
m3 <- apply(XXXX^3 * yyyy,2,function(x){sum(x)*step})
m4 <- apply(XXXX^4 * yyyy,2,function(x){sum(x)*step})
print(m4)

T <- dim(y)[1]
n <- dim(y)[2]

B_chol <- t(chol(Sigma_eta)) 
u <- Eta %*% t(solve(B_chol)) 
Y <- u  #get unit variance residuals

# RUN ICA----------------------------------------------------
best.value <- 100000000000000000
all.permut <- do.permut(n)
for(permut in 1){
  #for(permut in 1:dim(all.permut)[3]){
  # Run numerical optimizations:
  AA.permut <- solve(all.permut[,,permut]*1.01 + diag(n)) %*% (all.permut[,,permut] - diag(n))
  AA.0 <- AA.permut[lower.tri(AA.permut)]
  res.optim <- optim(AA.0,func.2.minimize,
                     gr = d.func.2.minimize,
                     method="Nelder-Mead",
                     # method="BFGS",
                     # method="CG",
                     control=list(trace=FALSE,maxit=1000))
  AA.0 <- res.optim$par
  res.optim <- optim(AA.0,func.2.minimize,d.func.2.minimize,
                     # method="Nelder-Mead",
                     method="BFGS",
                     # method="CG",
                     control=list(trace=FALSE))
  AA.0 <- res.optim$par
  res.optim <- optim(AA.0,func.2.minimize,
                     gr = d.func.2.minimize,
                     method="Nelder-Mead",
                     # method="BFGS",
                     # method="CG",
                     control=list(trace=FALSE,maxit=1000))
  print(res.optim$value)
  if(res.optim$value < best.value){
    AA.best <- res.optim$par
    best.value <- res.optim$value
  }
}

AA.est <- AA.best
n <- ncol(Y)
M <- make.M(n)
A.est <- matrix(M %*% AA.est,n,n)
C.PML <- (diag(n) + A.est) %*% solve(diag(n) - A.est)
print(C.PML)
C.PML%*%t(C.PML)

eps.PML <- Y %*% C.PML

A <- make.A.matrix(eps.PML,distri,C.PML)
eigen(A)$value
Omega <- make.Omega(eps.PML,distri)
eigen(Omega)$value

# Compute asymptotic covariance matrix of C.PML:
V <- make.Asympt.Cov.delta(eps.PML,distri,C.PML)

# collect estimation results
res.estim <- 
  cbind(
    c(C.PML),
    sqrt(diag(V)),
    c(C.PML)/sqrt(diag(V))
)

# Compute IC-IS measure 
shares<-t(Psi1[1,])%*%C.PML
shares2<-shares^2
rw.var<-t(Psi1[1,])%*%(C.PML%*%t(C.PML))%*%Psi1[1,]
IC.IS<-c(shares2[1]/rw.var,shares2[2]/rw.var,shares2[3]/rw.var,shares2[4]/rw.var)
IC.IS

#IC.IS.excevents <- IC.IS


# Wald tests
i_seq <- seq(from=0, to=16, by=4)

S_i <- NULL
S_ivar <- NULL 
W_i <- NULL
pvalues <- NULL
for (i in 1:4) {
  S_i[i] <- (t(Psi1[1,])%*%C.PML[,i])^2
  S_ivar[i] <- sum(res.estim[,2][(i_seq[i]+1):i_seq[i+1]])
}
W_i <- S_i/S_ivar  # Wald statistics
pvalues <- pchisq(W_i, df=1, lower.tail = FALSE)
Wald_tests <- cbind(S_i, S_ivar, W_i, pvalues)

# Compute IC-IS distribution (Non-central Beta distribution)-------------

# Each IC-IS_j is a non-central Beta Dist.
# psi_i*c_i is a non-central chi-square distribution, with non central parameter lambda = Sum mu_i  
# I compute the non-central Betas being the ratio of non-central chi-square distributions
library(sadists) #package needed for the doubly non-central version of the Beta

N <- 10000 # n° of draws
lambda_i <- shares2  # non-central parameters
dof <- c(1,3) # degrees of freedom of the X^2

dbeta <- matrix(0, N, n)
for (i in 1:n) {
  dbeta[,i] <- rdnbeta(N, df1 = dof[1], df2 = dof[2], ncp1 =lambda_i[i], ncp2 = sum(lambda_i[-i]))
}

beta <- matrix(0, N, n)
for (i in 1:n) {
  beta[,i] <- rbeta(N, shape1 = 1/2, shape2 = 3/2, ncp = lambda_i[i])
}

palette <- c(1,2,4,8)
for (i in 1:n) {
  x <- beta[,i]
  if(i==1) {    
    plot(density(x,bw=0.05, from = 0),col=palette[i],ylim=c(0,5), xlim=c(.0,.95), main = "IC-IS estimates and distributions", cex.main=1)
    abline(v=IC.IS[i],col=palette[i], lwd=2.5)
    }
  else {
    lines(density(x,bw=0.05, from = 0),col=palette[i])
    abline(v=IC.IS[i],col=palette[i], lwd=2.5)
    }
}
legend("topright", legend=c("Participant-Bid", "Participant-Ask", "Sip-Bid", "Sip-Ask"),
       col=c(1,2,4,8), lty=1, lwd = 2.5, cex=0.8)


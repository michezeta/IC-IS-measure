#------ Set simulation parameters ------#

n <- 4     # cross-sectional dimension
T <- 5000  # time series dimension
M <- 500   # number of simulations

# diurnal U-shape pattern for the variance (these parameters are fixed)
C <- 1
A <- 0.75
B <- 0.25
a <- 10
b <- 10 
s_t <- NULL

t <- seq(from=0, to=1, length.out = T)
for (i in 1:length(t)) {
  s_t[i] <- C + A*exp(-a*t[i]) + B*exp(-b*(1-t[i]))
}
s_t2 <- s_t^2  # variance following the U-shape pattern
plot(s_t2, type="l")

# VECM coefficients (I use just one lag for simplicity , captured by phi1)

phi1=matrix(c(0.4,0.6,0.2,0.1,-0.9,0.35,-0.2,0.35,-0.25,0.55,-0.7,0.6,0.3,-0.1,0.4,0.1), nrow=4, ncol=4, byrow = FALSE)
alpha<-matrix(c(0.09,0.1,0.025,0.08,0.06,0.01,0.05,0.07,0.09,0.04,0.03,0.06), nrow=4,ncol=3)
Beta<-matrix(c(1,-1,0,0,1,0,-1,0,1,0,0,-1),nrow = 4,ncol = 3)

# Mixing matrices (Choleski to whiten the innovations and orthogonal C) 

C <- matrix(c(0.99578180,-0.05120370,-0.05020271,0.05656933,0.05554580,0.99659534,  
              0.05539517,-0.02685817,0.04791627,-0.05699772,0.99712765,
              -0.01064394,-0.05582357,0.02914462,0.01498762,0.99794149), ncol = 4, nrow = 4) # true orthogonal matrix C
S <- matrix(c(0.9,0.4,0.5,0.3,0,0.6,0.2,0.5,0,0,0.7,0.3,0,0,0,0.1), nrow = 4, ncol = 4) # true Choleski dec. matrix S


#=======================================#
# ---- Montecarlo simulation begins ----#
#=======================================#
library(mcompanion)

IC.IS_est <- matrix(0, nrow = M, ncol = n)

for (m in 1:M) {
  
# --- Generate shocks from Student distributions ---
  
  v_t <- (2*s_t2)/(s_t2 - 1)      #time-varying degree of freedom of the Student
  vs_t <- cbind(s_t2, v_t)
  
  epsilons <- matrix(0, nrow = T, ncol = n)
  eps <- NULL
  for (j in 1:n) {
    for (i in 1:T) {
      eps[i] <- rt(1, df=v_t[i])
    }
    epsilons[,j] <- eps
  }  

# ---------- generate VECM model ---------
  
  u <- t(S%*%C%*%t(epsilons)) # cross-correlated price innovations
  Omega <- cov(u) # covariance matrix of the correlated u_t
  
  # set initial values for the simulation
  Pt <- matrix(0, nrow = T, ncol = 4)
  dPt <- matrix(0, nrow = T, ncol = 4)
  dPt[1,] <- c(0.05, 0.07, 0.09,0.08) # initial values for dP0
  Pt[1,] <- c(10.5,10,9,9.5) # initial values for P0
  
  for (i in 1:(T-1)) {
    dPt[i+1,] <- alpha%*%t(Beta)%*%Pt[i,] + phi1%*%dPt[i,] + u[i+1,]
    Pt[i+1,] = Pt[i,] + dPt[i+1,]
  }
    
# ----- Estimate the VECM, get residuals ----
  
  nb.lags = 1  # always 1 lag for simplicity, do not set this parameter differently!
  B <- cbind(matrix(c(1,1,1), ncol = 1,  nrow = dim(Y)[2]-1), - diag(dim(Y)[2]-1)) # cointegrating matrix is known
  
  X <- NULL
  for(i in 1:nb.lags){
    lagged.Y <- rbind(
      matrix(NaN,i,n),
      Pt[1:(T-i),])
    X <- cbind(X,lagged.Y)
  }
  y <- diff(Pt, lag = 1) 
  x <- diff(X, lag = 1)
  BPt <- t(B%*%t(Pt))[-nrow(Pt),]
  
  Phi <- matrix(0,n,n*nb.lags)
  a <- matrix(0,n,n-1)
  Eta <- NULL # Eta is the matrix of OLS residuals
  
  # estimate equations and get residuals
  for(i in 1:n){
    eq <- lm(y[,i] ~ BPt + x + 0)
    Phi[i,] <- eq$coef[n:(dim(Phi)[2]+(n-1))] 
    a[i,] <- eq$coef[1:(n-1)]
    Eta <- cbind(Eta,eq$residuals)
  }

# --- Apply the PML procedure to estimate the mixing matrix C ---  
  
  distri <- list(
    type=c("student","student","student", "student"),
    df=c(4,5,7,12),
    p=c(.5,.5,.5,.5),
    mu=c(.1,.1,.1,.1),
    sigma=c(.5,.7,.9,1.1)
  )
  
  B_chol <- t(chol(cov(Eta))) 
  Y <- Eta %*% t(solve(B_chol)) 
  
# ------------ RUN ICA -----------
  
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
    #print(res.optim$value)
    if(res.optim$value < best.value){
      AA.best <- res.optim$par
      best.value <- res.optim$value
    }
  }
  
  AA.est <- AA.best
  n <- ncol(Y)
  M <- make.M(n)
  A.est <- matrix(M %*% AA.est,n,n)
  C.PML <- (diag(n) + A.est) %*% solve(diag(n) - A.est) # This is the matrix of interest for the IC-IS
  #print(C.PML)
  #C.PML%*%t(C.PML)
  eps.PML <- Y %*% C.PML
  
# ---- ESTIMATE IC-IS --- #
  
  a_orth <- null_complement(a)
  B_orth <- null_complement(t(B))
  I <- diag(n)
  H_est <- solve(t(a_orth)%*%(I-Phi)%*%B_orth)
  Psi1_est <- B_orth%*%H_est%*%t(a_orth)
  
  shares_est<-t(Psi1_est[1,])%*%C.PML
  shares2_est<-shares_est^2
  rw.var_est<-as.numeric(t(Psi1_est[1,])%*%(C.PML%*%t(C.PML))%*%Psi1_est[1,])
  IC.IS_est[m,]<-shares2_est/rw.var_est
}

ICIS.mean <- c(mean(IC.IS_est[,1]),mean(IC.IS_est[,2]),mean(IC.IS_est[,3]),mean(IC.IS_est[,4]))
ICIS.sd <- c(sd(IC.IS_est[,1]),sd(IC.IS_est[,2]),sd(IC.IS_est[,3]),sd(IC.IS_est[,4]))
ICIS_results <- rbind(ICIS.mean,ICIS.sd)


# ---------------------- end of Montecarlo analysis ---------------------------------
# now I compare the estimated ICIS with the true ones and with the Choleski approach
#====================================================================================

#--- TRUE IC-IS ---#
#library(mcompanion)

alpha_orth <- null_complement(alpha)
B_orth <- null_complement(t(B))
I <- diag(n)
H <- solve(t(alpha_orth)%*%(I-phi1)%*%B_orth)
Psi1 <- B_orth%*%H%*%t(alpha_orth)

shares<-t(Psi1[1,])%*%C
shares2<-shares^2
rw.var<-as.numeric(t(Psi1[1,])%*%(C%*%t(C))%*%Psi1[1,])
IC.IS<-shares2/rw.var
IC.IS

#---- plot simulation results ----
palette <- c(1,2,4,8)

plot(density(IC.IS_est[,4], bw=0.05),xlim=c(.0,.9),ylim=c(.0,8), main = "(d)", xlab = "IC-IS", col=palette[4])
abline(v=IC.IS[,4], col=palette[4], lty= 2, lwd=2)
abline(v=ICIS.mean[4], col=palette[4])
for (i in 1:3) {
  lines(density(IC.IS_est[,i], bw=0.05), col=palette[i])
  abline(v=IC.IS[,i], col=palette[i], lty= 2, lwd=2)
  abline(v=ICIS.mean[i], col=palette[i])
}
legend("topright", legend=c("estimated IC-IS", "True IS"),
       col=c(1,1), lty=c(1,2), lwd = 2.5, cex=0.8)



#---------  all Choleski permutation based IS ----------------#
# no Montecarlo simulation here for the moment, not necessary #

library(combinat)
perms <- permn(c(1,2,3,4))

IC.IS_mat <- matrix(0,nrow = length(perms), ncol = n)
for (p in 1:length(perms)) {
  perm.eta <- Eta[,c(perms[[p]])]  # here chose the permutation of the variables in Choleski
  C_chol <- t(chol(cov(perm.eta))) # 
  Psi1_perm <- Psi1_est[,c(perms[[p]])] # permute columns of Psi1 accordingly
  rw.var_est <-as.numeric(t(Psi1_perm[1,])%*%(cov(perm.eta))%*%Psi1_perm[1,])
  shares_perm <- t(Psi1_perm[1,])%*%C_chol
  shares2_perm <- shares_perm^2
  IC.IS_perm <- rbind(shares2_perm/rw.var_est, perms[[p]])
  IC.IS_perm <- IC.IS_perm[,order(IC.IS_perm[2,], decreasing = FALSE)]
  IC.IS_mat[p,] <- IC.IS_perm[1,]
}
IC.IS_min <- c(min(IC.IS_mat[,1]),min(IC.IS_mat[,2]),min(IC.IS_mat[,3]),min(IC.IS_mat[,4]))
IC.IS_max <- c(max(IC.IS_mat[,1]),max(IC.IS_mat[,2]),max(IC.IS_mat[,3]),max(IC.IS_mat[,4]))
IC.IS_minmax <- rbind(IC.IS_min, IC.IS_max)
IC.IS_minmaxmid <- rbind(IC.IS_minmax,colMeans(IC.IS_minmax[,c(1,2,3,4)]))






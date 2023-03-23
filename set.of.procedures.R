



simul.distri <- function(distri,nb.sim){
  # Simulation of independent shocks
  eps <- NULL
  for(i in 1:length(distri$type)){
    if(distri$type[i]=="gaussian"){
      eps.i <- rnorm(nb.sim)
    }else if(distri$type[i]=="mixt.gaussian"){
      p <- distri$p[i]
      mu.1 <- distri$mu[i]
      sigma.1 <- distri$sigma[i]
      mu.2 <- - p/(1-p)*mu.1
      sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )
      B <- (runif(nb.sim)<p)*1
      eps.i <- B*rnorm(nb.sim,mean = mu.1,sd = sigma.1) +
        (1-B)*rnorm(nb.sim,mean = mu.2,sd = sigma.2)
    }else if(distri$type[i]=="student"){
      nu <- distri$df[i]
      eps.i <- rt(nb.sim,df = nu)/sqrt(nu/(nu-2))
    }else if(distri$type[i]=="laplace"){
      U <- runif(nb.sim) - .5
      b <- 1/sqrt(2)
      eps.i <- - b * sign(U) * log(1 - 2 * abs(U))
    }else if(distri$type[i]=="hyper.sec"){
      U <- runif(nb.sim)
      eps.i <- 2/pi * log(tan(pi/2*U))
    }
    eps <- cbind(eps,eps.i)
  }
  return(eps)
}



log.g.gaussian <- function(x,mu,sigma,indic.deriv=0){
  # Gaussian distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  log.g <- -(x - mu)^2/(2*sigma^2)
  if(indic.deriv==1){
    d.log.g <- -(x - mu)/sigma^2
    d2.log.g <- -1/sigma^2 * rep(1,length(x))
  }else{
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g=log.g,d.log.g=d.log.g,d2.log.g=d2.log.g))
}


log.g.student <- function(x,nu,indic.deriv=0){
  # Student distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  log.g <- -(1 + nu)/2 * log(1 + (x * sqrt(nu/(nu-2)))^2/nu)
  if(indic.deriv==1){
    d.log.g <- - x * (1 + nu) / (nu - 2 + x^2)
    d2.log.g <- - (1 + nu) * (nu - 2 - x^2) / (nu - 2 + x^2)^2
  }else{
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g=log.g,d.log.g=d.log.g,d2.log.g=d2.log.g))
}

log.g.laplace <- function(x,indic.deriv=0){
  # Student distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  b <- 1/sqrt(2)
  log.g <- - abs(x)/b
  if(indic.deriv==1){
    d.log.g <- - sign(x) / b
    d2.log.g <- 0
  }else{
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g=log.g,d.log.g=d.log.g,d2.log.g=d2.log.g))
}


log.g.hyper.sec <- function(x,indic.deriv=0){
  # Student distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  log.g <- - log(exp(pi/2*x)+exp(-pi/2*x))
  if(indic.deriv==1){
    d.log.g <- - (1/(exp(pi/2*x)+exp(-pi/2*x))*(exp(pi/2*x)*pi/2+exp(-pi/2*x)*(-pi/2)))
    d2.log.g <- ((-pi^2*exp(-1/2*pi*x)^2)/(exp(-1/2*pi*x)^4+2*exp(-1/2*pi*x)^2+1))
  }else{
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g=log.g,d.log.g=d.log.g,d2.log.g=d2.log.g))
}


func.2.min.rec <- function(c1,Y,distri){
  # This function is used in the context of recursive PML
  c1.used <- sign(c1)*pmin(abs(c1),1)
  c2 <- sqrt(1 - c1.used^2)
  c.i <- rbind(c1.used,c2)
  eps <- Y %*% c.i
  if(distri$type == "student"){
    g.aux <- log.g.student(eps,nu = distri$df)
  }else if(distri$type == "gaussian"){
    g.aux <- log.g.gaussian(eps,mu = 0,sigma = 1)
  }else if(distri$type == "mixt.gaussian"){
    p <- distri$p
    mu.1 <- distri$mu
    sigma.1 <- distri$sigma
    mu.2 <- - p/(1-p)*mu.1
    sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )
    g.aux <- log.g.mixt.gaussian(eps,mu = c(mu.1,mu.2),sigma = c(sigma.1,sigma.2),p)
  }else if(distri$type == "laplace"){
    g.aux <- log.g.laplace(eps)
  }else if(distri$type == "hyper.sec"){
    g.aux <- log.g.hyper.sec(eps)
  }
  penalty <- 1000000*(abs(c1)-1)*(abs(c1)>1)
  return(-apply(g.aux$log.g,2,sum) + penalty)
}


make.M <- function(n){
  # function that builds M_n, the n2 x ((n-1)n/2) matrix that is such that Vec(A) = M_n a where a is the vectorized entries of the
  # lower triangular matrix of A and A is such that t(A) = -A.
  aux1 <- 1:(n^2)
  aux2 <- matrix(0,n,n)
  aux2[lower.tri(aux2)] <- 1:((n-1)*n/2)
  aux2 <- matrix(aux2,ncol=1)
  aux3 <- matrix(0,n,n)
  aux3[lower.tri(aux3)] <- 1:((n-1)*n/2)
  aux3 <- t(aux3)
  aux3 <- matrix(aux3,ncol=1)
  
  M <- matrix(0,n^2,n*(n-1)/2)
  indic.row <- which(aux2>0)
  M[indic.row,] <- diag(n*(n-1)/2)
  indic.row <- which(aux3>0)
  M[indic.row,] <- - diag(n*(n-1)/2)
  return(M)
}


make.C <- function(AA,n){
  A <- matrix(0,n,n)
  A[lower.tri(A)] <- AA
  A[t(lower.tri(A))] <- - AA
  if(prod(A==(0*A))==1){
    C <- diag(n)
  }else{
    C <- (diag(n) + A) %*% solve(diag(n) - A)
  }
  return(C)
}

make.AA <- function(C){
  n <- dim(C)[1]
  A <- matrix(0,n,n)
  Id <- diag(n)
  A <- solve(C + Id) %*% (C - Id)
  return(A)
}

pseudo.log.L <- function(Y,AA,distri,indic.Jacobian=0){
  # AA is \mathcal{A}
  # Y is the matrix of observations
  # g is a likelihood function
  n <- ncol(Y)
  M <- make.M(n)
  I <- diag(n)
  A <- matrix(M %*% AA,n,n)
  if (n==4){
    A[2,3] <- - A[3,2]
    A[1,4] <- - A[4,1]
  }
  C <- make.C(AA,n)
  eps <- Y %*% C
  g.eps <- g(eps,distri,indic.deriv = indic.Jacobian)
  vector.log.L <- g.eps$log.L
  log.L <- sum(vector.log.L)
  if(indic.Jacobian==0){
    Jacobian <- NULL
  }else{
    aux1 <- (t(solve(I - A)) %x% (I + A)) %*% M
    d.log.L <- g.eps$d.log.L # This is a T x n matrix
    aux2.1 <- matrix(1,1,n) %x% Y
    aux2.2 <- d.log.L %x% matrix(1,1,n)
    aux2 <- aux2.1 * aux2.2
    Jacobian <- apply(aux2 %*% aux1,2,sum)
  }
  return(list(log.L=log.L,Jacobian=Jacobian))
}

pseudo.log.L.used <- function(AA,Y,distri){
  res <- pseudo.log.L(Y,AA,distri)
  return(-res$log.L)
}

g <- function(eps,distri,indic.deriv=0){
  # Y is of dimension T x n
  # distri is a list:
  #       - distri$type contains the types of distributions ("student" or "gaussian")
  #       - distri$df contains the degrees of freedom for student distri (NaN for non student distri)
  n <- dim(eps)[2]
  log.L <- NULL
  d.log.L <- NULL
  d2.log.L <- NULL
  for(ii in 1:n){
    if(distri$type[ii] == "student"){
      g.aux <- log.g.student(eps[,ii],nu = distri$df[ii],indic.deriv)
    }else if(distri$type[ii] == "gaussian"){
      g.aux <- log.g.gaussian(eps[,ii],mu = 0,sigma = 1,indic.deriv)
    }else if(distri$type[ii] == "mixt.gaussian"){
      p = distri$p[ii]
      mu.1 <- distri$mu[ii]
      sigma.1 <- distri$sigma[ii]
      mu.2 <- - p/(1-p)*mu.1
      sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )
      g.aux <- log.g.mixt.gaussian(eps[,ii],mu = c(mu.1,mu.2),sigma = c(sigma.1,sigma.2),
                                   p, indic.deriv)
    }else if(distri$type[ii] == "laplace"){
      g.aux <- log.g.laplace(eps[,ii], indic.deriv)
    }else if(distri$type[ii] == "hyper.sec"){
      g.aux <- log.g.hyper.sec(eps[,ii], indic.deriv)
    }
    log.L    <- cbind(log.L,g.aux$log.g)
    d.log.L  <- cbind(d.log.L,g.aux$d.log.g)
    d2.log.L <- cbind(d2.log.L,g.aux$d2.log.g)
  }
  return(list(log.L=log.L,d.log.L=d.log.L,d2.log.L=d2.log.L))
}

func.2.minimize <- function(theta){
  aux <- pseudo.log.L(Y,theta,distri)
  return(- aux$log.L)
}

d.func.2.minimize <- function(theta){
  aux <- pseudo.log.L(Y,theta,distri,indic.Jacobian = 1)
  return(- aux$Jacobian)
}

make.g.stars <- function(eps,distri){
  res.g <- g(eps,distri,indic.deriv=1)
  g.1.star <- apply(res.g$d.log.L,2,mean)
  g.2.star <- apply(res.g$d.log.L^2,2,mean)
  g.3.star <- apply(eps*res.g$d.log.L,2,mean)
  g.4.star <- apply(res.g$d2.log.L,2,mean)
  return(list(g.1.star=g.1.star,g.2.star=g.2.star,g.3.star=g.3.star,g.4.star=g.4.star))
}

make.Omega <- function(eps,distri){
  n <- dim(eps)[2]
  
  all.g.stars <- make.g.stars(eps,distri)
  g.1 <- matrix(all.g.stars$g.1.star,ncol=1)
  g.2 <- matrix(all.g.stars$g.2.star,ncol=1)
  g.3 <- matrix(all.g.stars$g.3.star,ncol=1)
  vec.1 <- matrix(1,n,1)
  
  g.1.g.1 <- g.1 %*% t(g.1)
  
  Omega <- g.1.g.1 %x% diag(n)
  
  for(i in 1:(n-1)){
    aux <- diag(n) %x% g.1.g.1[,i]
    Omega[(i*n+1):n^2,((i-1)*n+1):(i*n)] <- 
      Omega[(i*n+1):n^2,((i-1)*n+1):(i*n)] - aux[(i*n+1):n^2,]
  }
  
  # Add lower-tri parts of block diagonal:
  Omega <- Omega + diag(n) %x% g.1.g.1
  
  # Remove all entries on and above diag:
  diag(Omega) <- 0
  Omega[upper.tri(Omega)] <- 0
  
  # Make it symmetric:
  Omega <- Omega + t(Omega)
  
  # Fill diagonal:
  aux <- g.2 %*% t(vec.1) + vec.1 %*% t(g.2) - 2 * (g.3 %*% t(vec.1)) * (vec.1 %*% t(g.3))
  diag(Omega) <- c(aux)
  
  # Remove rows and columns to get a [n(n-1)/2] x [n(n-1)/2] matrix:
  aux <- upper.tri(diag(n))
  diag(aux) <- TRUE
  indic.rows.to.keep <- which(!aux)
  
  Omega <- Omega[indic.rows.to.keep,]
  if(class(Omega)[1]=="matrix"){
    Omega <- Omega[,indic.rows.to.keep]    
  }else{
    Omega <- Omega[indic.rows.to.keep]
  }

  return(Omega)
}

make.A.matrix <- function(eps,distri,C){
  n <- dim(eps)[2]
  
  all.g.stars <- make.g.stars(eps,distri)
  g.1 <- matrix(all.g.stars$g.1.star,ncol=1)
  g.2 <- matrix(all.g.stars$g.2.star,ncol=1)
  g.3 <- matrix(all.g.stars$g.3.star,ncol=1)
  g.4 <- matrix(all.g.stars$g.4.star,ncol=1)
  vec.1 <- matrix(1,n,1)
  
  A.1 <- matrix(0,n^2,n^2)
  A.2 <- matrix(0,n^2,n^2)
  for(i in 1:n){
    # compute the a_{i,.}:
    A_i <- (-g.4[i] + g.3 %*% t(vec.1)) * t(C)
    # Put these in A:
    A.1[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)] <- A_i
    A.2[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)] <- t(C)
  }
  for(j in 1:n){
    # compute the a_{.,j}:
    A_j <- (-g.4 %*% t(vec.1) + g.3[j]) * (vec.1 %*% matrix(C[,j],nrow=1))
    # Put these in A (diagonal blocks):
    for(i in 1:n){
      A.1[((j-1)*n+i),((i-1)*n+1):(i*n)] <- - A_j[i,]
      A.2[((j-1)*n+i),((i-1)*n+1):(i*n)] <- c(C[,j])
    }
  }
  # remove rows from A.1
  aux <- upper.tri(diag(n))
  diag(aux) <- TRUE
  indic.rows.to.keep <- which(!aux)
  A.1 <- A.1[indic.rows.to.keep,]
  A.2 <- A.2[indic.rows.to.keep,]
  
  A.3 <- matrix(0,n,n^2)
  for(i in 1:n){
    A.3[i,((i-1)*n+1):(i*n)] <- c(C[,i])
  }
  
  A <- rbind(A.1,A.2,A.3)
  
  return(A)
}


make.Asympt.Cov.delta <- function(eps,distri,C){
  n <- dim(eps)[2]
  T <- dim(eps)[1]
  A <- make.A.matrix(eps,distri,C)
  Omega <- make.Omega(eps,distri)
  Omega.augment <- matrix(0,n^2,n^2)
  Omega.augment[1:(n*(n-1)/2),1:(n*(n-1)/2)] <- Omega
  A_1 <- solve(A)
  mat.cov <- 1/T * A_1 %*% Omega.augment %*% t(A_1)
  return(mat.cov)
}

make.Asympt.Cov.deltaAronde <- function(eps,distri,C){
  A0 <- make.AA(C)
  n <- dim(eps)[2]
  T <- dim(eps)[1]
  Id <- diag(n)
  M.n <- matrix(0,n*(n-1)/2,n^2)
  aux <- matrix(1:(n^2),n,n)
  indic.lower.tri <- aux[lower.tri(aux)]
  for(i in 1:length(indic.lower.tri)){
    M.n[i,indic.lower.tri[i]] <- 1
  }
  mat.cov.C <- make.Asympt.Cov.delta(eps,distri,C)
  kron.mat <- t(Id - A0) %x% (Id - A0)
  mat.cov.Aronde <- 1/4 * M.n %*% kron.mat %*% mat.cov.C %*% t(kron.mat) %*% t(M.n)
  return(mat.cov.Aronde)
}

make.Asympt.Cov.delta_1 <- function(eps,distri,C){
  n <- dim(eps)[2]
  T <- dim(eps)[1]
  A <- make.A.matrix(eps,distri,C)
  Omega <- make.Omega(eps,distri)
  Omega.augment <- matrix(0,n^2,n^2)
  Omega.augment[1:(n*(n-1)/2),1:(n*(n-1)/2)] <- solve(Omega)
  mat.cov <- T * t(A) %*% Omega.augment %*% A
  return(mat.cov)
}



do.permut <- function(n){
  if(n==1){
    all.permut <- array(1,c(1,1,1))
  }else{
    permut_n_1 <- do.permut(n-1)
    all.permut <- array(NaN,c(n,n,factorial(n)))
    for(i in 1:n){
      aux <- array(NaN,c(n,n,factorial(n-1)))
      aux[,1,] <- 0
      aux[i,,] <- 0
      aux[i,1,] <- 1
      if(i>1){
        aux[(1:(i-1)),2:n,] <-  permut_n_1[(1:(i-1)),,]
      }
      if(i<n){
        aux[(i+1):n,2:n,] <-  permut_n_1[i:(n-1),,]
      }
      all.permut[,,((i-1)*factorial(n-1)+1):(i*factorial(n-1))] <- aux
    }
  }
  return(all.permut)
}


do.signs <- function(A){
  n <- dim(A)[2]
  #print(n)
  if(n==1){
    A.all.sign <- array(NaN,c(dim(A)[1],1,2*dim(A)[3]))
    A.all.sign[,,1:dim(A)[3]] <- A
    A.all.sign[,,(dim(A)[3]+1):(2*dim(A)[3])] <- -A
  }else{
    A.all.sign <- array(NaN,c(dim(A)[1],dim(A)[2],2^n*dim(A)[3]))
    aux <- array(A[,2:dim(A)[2],],c(dim(A)[1],dim(A)[2]-1,dim(A)[3]))
    A.n_1 <- do.signs(aux)
    aux <- array(NaN,c(dim(A)[1],dim(A)[2],2^(n-1)*dim(A)[3]))
    aux[,1,] <- A[,1,]
    aux[,2:dim(A)[2],] <- A.n_1
    A.all.sign[,,1:dim(A.n_1)[3]] <- aux
    aux <- array(NaN,c(dim(A)[1],dim(A)[2],2^(n-1)*dim(A)[3]))
    aux[,1,] <- -A[,1,]
    aux[,2:dim(A)[2],] <- A.n_1
    A.all.sign[,,(dim(A.n_1)[3]+1):(2*dim(A.n_1)[3])] <- aux
  }
  return(A.all.sign)
}


find.permut <- function(C,C.permut,distance="eucl"){
  n <- dim(C)[1]
  best.ddiff <- 10000000000
  all.permut <- do.permut(n)
  for(i in 1:dim(all.permut)[3]){
    C.aux <- C.permut %*% all.permut[,,i]
    C.aux <- array(C.aux,c(dim(C.aux)[1],dim(C.aux)[2],1))
    C.aux.diff.signs <- do.signs(C.aux)
    for(j in 1:dim(C.aux.diff.signs)[3]){
      if(distance=="eucl"){
        ddiff <- sum((C-C.aux.diff.signs[,,j])^2)
      }else{
        ddiff <- sum(abs(C-C.aux.diff.signs[,,j]))
      }
      if(ddiff < best.ddiff){
        C.permut.permut <- C.aux.diff.signs[,,j]
        best.ddiff <- ddiff
      }
    }
  }
  return(C.permut.permut)
}

make.nice.distri.name <- function(distri){
  vec.names <- NULL
  for(i in 1:length(distri$type)){
    if(distri$type[i]=="gaussian"){
      vec.names = c(vec.names,
                    "Gaussian")
    }else if(distri$type[i]=="student"){
      vec.names = c(vec.names,
                    paste("t(",toString(distri$df[i]),")",sep=""))
    }else if(distri$type[i]=="laplace"){
      vec.names = c(vec.names,"laplace")
    }else if(distri$type[i]=="hyper.sec"){
      vec.names = c(vec.names,"Hyperb. sec.")
    }else if(distri$type[i]=="mixt.gaussian"){
      vec.names = c(vec.names,"mixt.-Gaussian")
    }
  }
  return(vec.names)
}


round.fixed.length <- function(X,n){
  # This procedure is used in the automatic creation of Latex tables
  # x is numeric. The output is a string with n numbers after ".", even if they are 0s.
  signX <- sign(X)
  if(signX>0){
    signX <- NULL
  }else{
    signX <- "-"
  }
  x <- abs(X)
  string.integer <- toString(as.integer(x))
  string.decimal <- toString(round(x - as.integer(x,0),n))
  while(nchar(string.decimal)<n+2){
    string.decimal <- paste(string.decimal,"0",sep="")
  }
  if(n==0){
    string.x <- paste(signX,string.integer,sep="")
  }else{
    string.x <- paste(signX,string.integer,".",str_replace(string.decimal,"0.",""),sep="")
  }
  return(string.x)
}

log.g.mixt.gaussian <- function(x,mu,sigma,p,indic.deriv=0){
  # Mixture of Gaussian distributions
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  # mu is 2 x 1 vector; sigma is a 2 x1 vector, p is a scalar
  mu_1 <- mu[1]
  mu_2 <- mu[2]
  sigma_1 <- sigma[1]
  sigma_2 <- sigma[2]
  log.g <- log(
    p / sqrt(2*pi) / sigma_1 * exp(- (x - mu_1)^2/(2 * sigma_1^2)) +
      (1-p) / sqrt(2*pi) / sigma_2 * exp(- (x - mu_2)^2/(2 * sigma_2^2))
  )
  if(indic.deriv==1){
    # computation done on Xcas (online)
    d.log.g <- ((-mu_1*p*sigma_2^3*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))+
                   mu_2*p*sigma_1^3*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                   mu_2*sigma_1^3*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                   p*sigma_1^3*x*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+
                   p*sigma_2^3*x*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))+
                   sigma_1^3*x*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2)))
                /(p*sigma_1^3*sigma_2^2*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))
                  -p*sigma_1^2*sigma_2^3*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))
                  -sigma_1^3*sigma_2^2*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))))
    d2.log.g <- ((-mu_1^2*p^2*sigma_2^4*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*
                    exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+mu_1^2*p*sigma_2^4*
                    exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*
                    exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+2*mu_1*mu_2*p^2*sigma_1^2*sigma_2^2*
                    exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                    2*mu_1*mu_2*p*sigma_1^2*sigma_2^2*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*
                    exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-2*mu_1*p^2*sigma_1^2*sigma_2^2*x*
                    exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+
                    2*mu_1*p^2*sigma_2^4*x*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+
                    2*mu_1*p*sigma_1^2*sigma_2^2*x*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*
                    exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-2*mu_1*p*sigma_2^4*x*
                    exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                    mu_2^2*p^2*sigma_1^4*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+
                    mu_2^2*p*sigma_1^4*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+
                    2*mu_2*p^2*sigma_1^4*x*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                    2*mu_2*p^2*sigma_1^2*sigma_2^2*x*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                    2*mu_2*p*sigma_1^4*x*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+2*mu_2*p*sigma_1^2*sigma_2^2*
                    x*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-p^2*sigma_1^5*sigma_2*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))^2+
                    p^2*sigma_1^4*sigma_2^2*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-p^2*sigma_1^4*x^2*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*
                    exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+p^2*sigma_1^2*sigma_2^4*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+2*p^2*sigma_1^2*
                    sigma_2^2*x^2*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                    p^2*sigma_1*sigma_2^5*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))^2-p^2*sigma_2^4*x^2*
                    exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+2*p*sigma_1^5
                  *sigma_2*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))^2-p*sigma_1^4*sigma_2^2*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*
                    exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+p*sigma_1^4*x^2*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-
                    p*sigma_1^2*sigma_2^4*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-2*p*sigma_1^2*sigma_2^2*x^2*
                    exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+p*sigma_2^4*x^2*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*
                    exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))-sigma_1^5*sigma_2*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))^2)
                 /(p^2*sigma_1^5*sigma_2^3*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))^2-2*p^2*
                     sigma_1^4*sigma_2^4*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+p^2*sigma_1^3*sigma_2^5*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))^2-2*p*sigma_1^5*sigma_2^3*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))^2+2*p*sigma_1^4*
                     sigma_2^4*exp((-mu_1^2+2*mu_1*x-x^2)/(2*sigma_1^2))*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))+sigma_1^5*sigma_2^3*exp((-mu_2^2+2*mu_2*x-x^2)/(2*sigma_2^2))^2))
  }else{
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g=log.g,d.log.g=d.log.g,d2.log.g=d2.log.g))
}


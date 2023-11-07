update_random_effect <- function(flash){
  flash <- get.fit(flash)
  if(!uses.re(flash)) return(flash)

  cov.is.list <- inherits(flash[["re_cov"]], "list")

  n <- get.R2.n(flash)
  re_dim <- get.re.dim(flash)

  equal.D <- (n == 2 & re_dim == 1) | (re_dim == 2 & n == 1) | is.tau.constant(flash)

  ## calculate residuals
  EF <- get.EF(flash)
  ER <- get.Y(flash) - lowrank.expand(EF)
  if(re_dim == 2) ER <- t(ER)

  if(!cov.is.list & equal.D){
    D <- get.tau(flash)
    if(is.tau.constant(flash)) D <- rep(D, nrow(ER))
    re_cov <- get.re.cov(flash)
    ub <- update.b.all(ER, re_cov, D)
    flash$EB <- ub$EB
    flash$EB2 <- ub$EB2
    #flash <- set.KL(flash, ub$kl, n = 3)
    flash$KLB <- ub$kl
  }else{
    D <- get.tau(flash)
    if(is.tau.constant(flash)) D <- rep(D, ncol(ER))
    P <- nrow(ER)
    N <- ncol(ER)

    ebeb2 <- sapply(1:P, function(j){
      if(!equal.D) myD <- D[j]
        else myD <- D
      re_cov <- get.re.cov(flash, n = j)
      unlist(update.b.1(t(ER[j,,drop = FALSE]), re_cov, myD))
    }) |> matrix(nrow = P, byrow = TRUE)
    flash$EB <- ebeb2[,1:N]
    flash$EB2 <- ebeb2[, N+ (1:N)]
    #flash <- set.KL(flash, sum(ebeb2[,2*N + 1]), n = 3)
    flash$KLB <-  sum(ebeb2[,2*N + 1])
  }
  if(re_dim == 2){
    flash$EB <- t(flash$EB)
    flash$EB2 <- t(flash$EB2)
  }
  return(flash)

}

uses.re <- function(flash){
  if(is.null(flash[["re_cov"]])){
    return(FALSE)
  }
  return(TRUE)
}

get.re.cov <- function(f, n = NULL){
  re_cov <- f[["re_cov"]]
  if(is.null(n)) return(re_cov)
  if(inherits(re_cov, "list")){
    return(re_cov[[n]])
  }
  return(re_cov)
}


get.re.dim <- function(f){
  f[["re_dim"]]
}

get.EB <- function(f){
  EB <- f[["EB"]]
  if(is.null(EB)){
    d <- dim(get.Y(f))
    EB <- matrix(0, nrow = d[1], ncol = d[2])
  }
  return(EB)
}

get.EB2 <- function(f){
  EB2 <- f[["EB2"]]
  if(is.null(EB2)){
    d <- dim(get.Y(f))
    EB2 <- matrix(0, nrow = d[1], ncol = d[2])
  }
  return(EB2)
}

get.KLB <- function(f){
  KLB <- f[["KLB"]]
  if(is.null(KLB)) return(0)
  return(KLB)
}

## update b for one row/column
## eigS eigen(Sigma_i)
## D: diagonal elments of noise precision
## r residual Y - LF^T
## posterior:
# b | r \propto p(r | b)p(b)
# exp( (r - b)^T diag(D) (r-b) + b^T (S0^{-1}) b)
# b^T (S0^{-1} + diag(D)) b - 2*r^T diag(D) b
# Sb = solve((S0^{-1} + diag(D)) = U 1/(1/d0 + D) U^T
# ub = Sb diag(D) r
#
# KL calculation
# KL = E_post[ log(p_0(b)) - log(p_post(b))]
# p_0(b) = N(b; 0, S_0); S_0 = U diag(d_0) U^T; d_0 = eigS$values
# p_post(b) = N(b; mub, Sb); Sb = U diag(d) U^T; d = 1/(1/d0 + D)
# KL = log(1/sqrt(det(S_0))) - log(1/sqrt(det(Sb)))
#             - 0.5* E(b^T S0^{-1} b - (b - mub)^T Sb^{-1} (b - mub))
# E(b^T S0^{-1} b - (b - mub)^T Sb^{-1} (b - mub)) =
#      E[ b^T (S0^{-1} - Sb^{-1})b] + mub^T Sb^{-1} mub
#   = mub^T(S0^{-1} - Sb^{-1}) mub  + mub^T Sb^{-1} mub  + sum(diag(S0^{-1} - Sb^{-1})*diag(Sb))
#  =  mub^T U diag(1/d0) U^T mub +  sum(diag(S0^{-1} - Sb^{-1})*diag(Sb))

update.b.1 <- function(r, eigS, D){
  d0 <- eigS$values
  U <- eigS$vectors
  d <- 1/((1/d0) + D)
  # Sb <- emulator::quad.tform(diag(d), U) U diag(d) U^T
  Sb <- tcrossprod(U, tcrossprod(U, diag(d)))
  mub <- Sb %*% (r*D)
  mu2b <-  mub^2 + diag(Sb)
  ## calculate KL
  delta_logdet <- -0.5*sum(log(d0)) + 0.5*sum(log(d))

  # S0_inv <- emulator::quad.tform(diag(1/d0), U)
  S0_inv <- tcrossprod(U, tcrossprod(U, diag(1/d0)))
  # a <- emulator::quad.form(S0_inv , mub) mub^T S0_inv mub
  a <- crossprod(crossprod(S0_inv, mub), mub)
  X <- tcrossprod(U, tcrossprod(U, diag(-D)))
  b <- sum(diag(X)*diag(Sb))
  kl <- delta_logdet - 0.5*(a + b)

  return(list("EB" = mub, "EB2" = mu2b, "kl" = kl))
}

# for shared covariance matrix
update.b.all <- function(ER, eigS, D){
  d0 <- eigS$values
  U <- eigS$vectors
  d <- 1/((1/d0) + D)
  # Sb <- emulator::quad.tform(diag(d), U) U diag(d) U^T
  Sb <- tcrossprod(U, tcrossprod(U, diag(d)))
  EB <- t(Sb %*% (t(ER)*D))
  EB2 <- t(t(EB^2) + diag(Sb))

  ## calculate KL
  N <- nrow(ER)
  delta_logdet <- N*(-0.5*sum(log(d0)) + 0.5*sum(log(d)))

  # S0_inv <- emulator::quad.tform(diag(1/d0), U)
  S0_inv <- tcrossprod(U, tcrossprod(U, diag(1/d0)))
  # a <- sum(emulator::quad.tdiag(S0_inv , EB)) sum(mub^T S0_inv mub)
  a <- sum(rowSums(tcrossprod(EB, S0_inv) * EB))
  X <- tcrossprod(U, tcrossprod(U, diag(-D)))
  b <- N*(sum(diag(X)*diag(Sb)))
  kl <- delta_logdet - 0.5*(a + b)

  return(list(EB = EB, EB2 = EB2, kl = kl))
}

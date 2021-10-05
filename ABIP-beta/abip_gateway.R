require(MASS)
setwd("D:/Git/ABIP")
dyn.load("./ABIP-beta/src/abip.dll")

n <- 100
p <- 100
rho <- 0.5
Sigma_0 <- rho ^ abs(outer(1:p, 1:p, "-"))
Omega_0 <- solve(Sigma_0)
X <- mvrnorm(n, rep(0, p), Sigma_0)
Sigma_n <- t(X) %*% X / n

lambda <- c(1e-3, 1e-2, 1e-1)
nlambda <- length(lambda)
Omega <- rep(0, nlambda * p * p)
# x_list <- rep(0, 2*nlambda*p*p)

res <- .C("abip_R", as.double(Sigma_n), as.integer(p), as.double(lambda), as.integer(nlambda), Omega = as.double(Omega))
Omega <- matrix(res$Omega, p, nlambda * p)
print(Omega)

dyn.unload("./ABIP-beta/src/abip.dll")
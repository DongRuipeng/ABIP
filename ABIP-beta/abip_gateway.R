require(MASS)
setwd("D:/Git/ABIP")
dyn.load("./ABIP-beta/src/abip.dll")

n <- 100
p <- 4
rho <- 0.5
Omega_0 <- rho ^ abs(outer(1:p, 1:p, "-"))
Sigma_0 <- solve(Omega_0)
X <- mvrnorm(n, rep(0, p), Sigma_0)
Sigma_n <- t(X) %*% X / n 

lambda <- c(1,0.8)
nlambda <- 2

.C("abip_R", as.double(Sigma_n), as.integer(p), as.double(lambda), as.integer(nlambda))

dyn.unload("./ABIP-beta/src/abip.dll")
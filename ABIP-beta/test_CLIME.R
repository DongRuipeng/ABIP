require(MASS)
setwd("D:/Git/ABIP")

dll.path <- "./ABIP-beta/src/abip.dll"

n <- 100
p <- 100
rho <- 0.5
Sigma_0 <- rho ^ abs(outer(1:p, 1:p, "-"))
Omega_0 <- solve(Sigma_0)
X <- mvrnorm(n, rep(0, p), Sigma_0)
Sigma_n <- t(X) %*% X / n

lambda <- c(1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
lambda <- sort(lambda, decreasing = TRUE)
nlambda <- length(lambda)
Omega <- rep(0, nlambda * p * p)
time <- 0

dyn.load(dll.path)
res <- .C("abip_R", as.double(Sigma_n), as.integer(p), as.double(lambda), as.integer(nlambda), Omega = as.double(Omega), time = as.double(time))
dyn.unload(dll.path)

Omega <- matrix(res$Omega, p, nlambda * p)
time <- res$time

for (l in c(1:nlambda)) {
  error <- norm(Omega[, ((l - 1) * p + 1):(l * p)] - Omega_0, "F") / norm(Omega_0, "F")
  cat("The error of Omega[",l,"]: ",error,"\n", sep = "")
}
cat("the time is ",time," secs.\n")

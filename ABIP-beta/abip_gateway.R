dyn.load("./ABIP/abip_gateway.dll")
require(MASS)

p <- 200
rho <- 0.5
Sigma <- rho ^ abs(outer(1:p, 1:p, "-"))
result <- .C("sp_mat", as.double(Sigma), as.integer(p))

dyn.unload("./ABIP/abip_gateway.dll")
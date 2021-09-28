dyn.load("./Test/abip_r.dll")

p <- 3
rho <- 0.5 
Sigma <- rho^abs(outer(1:p,1:p,"-"))

result <- .C("print_mat", as.double(Sigma), as.integer(p))

dyn.unload("./Test/abip_r.dll")
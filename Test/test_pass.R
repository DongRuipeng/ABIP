dyn.load("./Test/test_pass.dll")

n <- 4
A <- matrix(c(1:n ^ 2), n, n)
B <- matrix(0, n, n)
C <- matrix(0, n, n)

result <- .C("add_mat", as.double(A), as.double(B), C = as.double(C), as.integer(n))

dyn.unload("./Test/test_pass.dll")
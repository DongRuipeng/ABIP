dyn.load("./Test/test_pass.dll")

n <- 4
A <- matrix(c(1:n ^ 2), n, n)
B <- A
C <- matrix(0, n, n)

result <- .C("add_mat", as.double(A), as.double(B), C = as.double(C), as.integer(n))

C <- matrix(result$C, n,n)
print(C)

dyn.unload("./Test/test_pass.dll")
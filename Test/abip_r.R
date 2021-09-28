dyn.load("./Test/abip_r.dll")

p <- 3
n <- p^2 
M <- matrix(c(1:n),p,p)
for(i in c(1:p)){
    for(j in c(1:p)){
        cat(M[i,j]," ")
    }
    cat("\n")
}
result <- .C("print_mat", as.double(M), as.integer(n))

dyn.unload("./Test/abip_r.dll")
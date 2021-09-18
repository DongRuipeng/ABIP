dyn.load("./Test/convolve.dll")

a <- rep(1, 5)
b <- rep(2, 5)

result <- .C("convolve",
a = as.double(a),
as.integer(length(a)),
b = as.double(b),
as.integer(length(b)),
ab = double(length(a) + length(b) - 1))
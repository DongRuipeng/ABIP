rho = 0.5;
p = 3;
Sigma = rho .^ abs((1:p) - (1:p)');
A = [Sigma, -Sigma, eye(p), zeros(p);...
    zeros(p), zeros(p), eye(p), eye(p)];
A = sparse(A);
test(A,27,13);
rho = 0.5;
n = 100;
p = 100;
Sigma_0 = rho .^ abs((1:p) - (1:p)');
Omega_0 = inv(Sigma_0);
X = mvnrnd(zeros(1,p),Sigma_0, n);
Sigma_hat = X' * X / n;

lambda = 1e-3;
Basics = eye(p);

likelihood = @(Sigma, Omega) ...
    (Sigma(:))' * Omega(:) - log(det(Omega));

A = [Sigma_hat, -Sigma_hat, eye(p), zeros(p);...
    zeros(p), zeros(p), eye(p), eye(p)];
A = sparse(A);
c = [ones(p,1); ones(p,1); zeros(p,1); zeros(p,1)];

data.A = A;
data.c = c;

% lbounds = zeros(4*p,1);
% ubounds = Inf(4*p,1);

% [A,b,c,info] = presolve(A,b,c,lbounds,ubounds);
% abips setting.
data.A = sparse(A);
data.c = full(c);
%     params_abips = struct('max_iters', 100000, 'max_outiters', 10000);
params_abips = struct("verbose",0);

% abips implementation.
tic;
for j=1:p
    b = [lambda * ones(p,1) + Basics(:,1);...
        2*lambda*ones(p,1)];
    data.b = full(b);
     [x, y, s, info_abips] = abip_direct(data, params_abips); 
end
time_abips = toc;
fprintf("the time is %f secs.\n", time_abips);
%     [x_abips, objp_abips] = postsolve(x, info);
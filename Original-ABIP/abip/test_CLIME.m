rho = 0.5;
n = 100;
p = 200;
Omega_0 = rho .^ abs((1:p) - (1:p)');
Sigma_0 = inv(Omega_0);
X = mvnrnd(zeros(1,p),Sigma_0, 2 * n);
X_train = X(1:n, :);
X_test = X((n+1): (2*n),:);
Sigma_hat = X_train' * X_train / n;

lambda = 1;
Basics = eye(p);

likelihood = @(Sigma, Omega) ...
    (Sigma(:))' * Omega(:) - log(det(Omega));

A = [Sigma_hat, -Sigma_hat, eye(p), zeros(p);...
    zeros(p), zeros(p), eye(p), eye(p)];
A = sparse(A);
b = [lambda * ones(p,1) + Basics(:,1);...
    2*lambda*ones(p,1)];
c = [ones(p,1); ones(p,1); zeros(p,1); zeros(p,1)];

data.A = A;
data.b = b;
data.c = c;

lbounds = zeros(4*p,1);
ubounds = Inf(4*p,1);

% [A,b,c,info] = presolve(A,b,c,lbounds,ubounds);

if ~info.feasible
    fprintf('The problem is infeasible!\n');
else
    % abips setting.
    data.A = sparse(A);
    data.c = full(c);
    data.b = full(b);
    params_abips = struct('max_iters', 100000, 'max_outiters', 10000);
    
    % abips implementation.
    tic; [x, y, s, info_abips] = abip_direct(data, params_abips); time_abips = toc;
%     [x_abips, objp_abips] = postsolve(x, info);
end
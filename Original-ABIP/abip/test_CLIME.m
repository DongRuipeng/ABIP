clc
rho = 0.5;
n = 100;
p = 80;
Sigma_0 = rho .^ abs((1:p) - (1:p)');
Omega_0 = inv(Sigma_0);
X = mvnrnd(zeros(1,p),Sigma_0, n);
Sigma_hat = X' * X / n;

lambda = [1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
nlambda = length(lambda);
Omega = zeros(p,p,nlambda);

A = [Sigma_hat, -Sigma_hat, eye(p), zeros(p);...
    zeros(p), zeros(p), eye(p), eye(p)];
A = sparse(A);
c = [ones(p,1); ones(p,1); zeros(p,1); zeros(p,1)];
Basics = eye(p);

data.A = A;
data.c = c;


% abips implementation.
tic;
for j=1:p
    params_abips = struct("verbose",0);
    for l = 1:nlambda
        b = [lambda(l) * ones(p,1) + Basics(:,j);...
            2*lambda(l)*ones(p,1)];
        data.b = full(b);
        [x, y, s, info_abips] = abip_direct(data, params_abips);
        params_abips = struct("verbose",0, "x",x,"y",y,"s",s);
        Omega(:,j,l) = x(1:p) - x((p+1):2*p);
    end
end
time_abips = toc;
for l=1:nlambda
    fprintf("The error of Omega[%i]: %f\n", ...
        l, norm(Omega(:,:,l)-Omega_0,'fro')...
        / norm(Omega_0,'fro'));
end
fprintf("the time is %f secs.\n", time_abips);

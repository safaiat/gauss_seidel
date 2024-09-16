close all;
clear, clc;

% generate A and b
resolution = 1;
[A,b] = generate_Ab(resolution);

% solving setting
maxit = 2000;
tol = 1;

% initialization
x = zeros(size(b));
r = zeros(size(x,1),maxit);

% Gaussian-Seidel iterations
tic;
for kk = 1:maxit
    x_k = x;
    for ii = 1:size(x,1)
        S1 = A(ii,1:ii-1)*x(1:ii-1);
        S2 = A(ii,ii+1:end)*x(ii+1:end);
        x(ii) = (b(ii)-S1-S2)/A(ii,ii);
    end
    x_k1 = x;
    r(:,kk) = b-A*x;

    if norm(r(:,kk),'fro') < tol
        iter_stop = kk;
        disp(['Converge at iteration ', num2str(kk)])
        break
    end

end
toc
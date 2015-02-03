
function [x,i,lambda, f, dualf, dualArg] = Newton2(A, b, c, T, x0, maxIter, threshold, alpha, beta)
%compute initial parameters
[m,n] = size(A);
x = zeros(n, maxIter);
x(:,1) = x0;

%begin iteration loop
for i = 1:maxIter
    %gradient and hessian
    f = T*c'*x(:,i) - sum(log(b-A*x(:,i)));
    d = 1./(b - A * x(:,i));
    grad = T*c + A' * d;
    Hess = A' * diag(d.^2) * A;
    
    %Newton step and decrement
    Dx = - Hess \ grad;
    lambda = - grad' * Dx;
    
    %stopping criterion
    if (lambda / 2) < threshold
        dualf = f -m/T;
        dualArg = 1/T * d;
        break
    else
        %backtracking line search
        t = 1;
        while T*c'*(x(:,i)+t*Dx) -sum(log(b-A*(x(:,i)+ t*Dx))) > (f+alpha*t*grad'*Dx)
            t = beta * t ;
        end 
        %update
        x(:,i+1) = x(:,i) + t * Dx ;
    end 
end 
end

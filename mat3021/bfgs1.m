%----------------------------------------------------------------------------------------------------
% An implementation of BFGS algorithm  for optimization problems
% Together with bracketing linesearch
%----------------------------------------------------------------------------------------------------

clear all;
close all;
format long

N=2; % number of variables (edit this according to objective function)
xi=sym(zeros(1,N));
for i=1:N
    syms("x"+i);
    xi(1,i) = ("x"+i);
end

% objective function
f = 100*((x2 - x1^2)^2) + (1 - x1)^2
x = [-1.2 1]; % Initial Guess
%x = [-0.5 0.2];
%f = 100*((x2 - x1^2)^2) + (1 - x1)^2 + 90*((x4 - x3^2)^2) + (1 - x3)^2 + 10.1*((x2 - 1)^2 + (x4 - 1)^2) - 19.8*(x2 - 1)*(x4 - 1)
%x = [-3 -1 -3 -1]; % Initial Guess

%f = (x1 + 10*x2)^2 + 5*(x3 - x4)^2 + (x2 - 2*x3)^4 + 10*(x1 - x4)^4
%x = [-3 -1 0 1]; % Initial Guess

% step size 
alpha=1;
c = 10^-4;
rho = 0.3;

I = eye(N);
e = 10^(-5); % Convergence Criteria
k = 1; % Iteration Counter

% Gradient Computation:
for i=1:N 
    df_dx(i)= diff(f, "x"+i);
end

G = subs(df_dx, xi, x);

% Hessian matrix Computation:
for i=1:N
    for j=1:N
      HF(i,j)=diff(diff(f,"x"+j), "x"+i); 
    end
end

HF_xk = subs(HF,xi,x);
%Pk = -1 * (HF_xk \ G');  % Search Direction
Hk = inv(HF_xk);
Pk = -Hk*G';
% display table
fprintf('k: %d\t', k);
for i=1:N
   fprintf('x%d: %d\t\t\t',i, x(i)); 
end
fprintf('f(xk): %d\t\t', subs(f,[xi],[x]));
fprintf('||∇f(xk)||: %d\t', norm(G)); 
fprintf('\n');

while norm(G) >= e
    x_old = x;
    % compute step length alpha_k
    
    lh = subs(f,xi, x + alpha * Pk');
    rh = subs(f,xi,x)+ c * alpha*G'*Pk';

    while lh>rh
        alpha = rho*alpha;
        lh = subs(f,xi,x+alpha*Pk');
        rh = subs(f,xi,x)+c*alpha*G'*Pk';
    end

    for i=1:N 
      x(i)= x(i) + alpha * Pk(i); % Updated xk value
    end
    % Update Hk 
    G_old = subs(df_dx, xi, x_old);
    G = subs(df_dx, xi, x);
    sk = x-x_old;
    yk = G - G_old;
    rhok = 1/ (sk*yk');
    Hk = (I-rhok * sk'* yk) * Hk * (I-rhok * yk' * sk) + rhok * sk' * sk;
    Pk = -1 * (Hk * G'); % New Search Direction
    k = k+1;
    % display table
    fprintf('k: %d\t', k); 
    for i=1:N
        fprintf('x%d: %d\t',i, x(i)); 
    end
    fprintf('f(xk): %d\t', double(subs(f,[xi], [x])));
    fprintf('||∇f(xk)||: %d\t', double(norm(G))); 
    fprintf('\n');

end
%----------------------------------------------------------------------------------------------------
% An implementation of Steepest descent method for optimization problems
% with bracketing line search method
%----------------------------------------------------------------------------------------------------

clear all;
close all;
format long

N=2;
xi=sym(zeros(1,N));
for i=1:N 
    syms("x"+i);
    xi(1,i) = ("x"+i);
    x(i) = 0; % Initial Guess
end
%x(1) =-0.5;
%x(2) = 0.2
x(1) =-1.2;
x(2) = 1.0;
f = 100*(x2-x1^2)^2+(1-x1)^2
%f = 5*x1^2 + x2^2 + 4*x1*x2 - 14*x1 - 6*x2 + 20
%f = 10*x1^2 + x2^2 + 5*x1*x2 - 14*x1 - 6*x2 + 10
%f = 10*x1^2 + 8*x2^2 + 5*x1*x2 - 14*x1 - 6*x2 + 10

e = 10^(-5); % Convergence Criteria
k = 1; % Iteration Counter
alpha= 1 ;
c = 10^-4;
rho = 0.3;

% Gradient Computation:
for i=1:N 
    df_dx(i)= diff(f, "x"+i);
end

G = subs(df_dx, xi, x);
Pk = -(G); % Search Direction

% display plot
figure(1); clf; fcontour(f, 'Fill', 'On'); axis equal; hold on
% display table
fprintf('k: %d\t', k);
for i=1:N
   fprintf('x%d: %d\t\t\t',i, x(i)); 
end
fprintf('f(xk): %d\t\t', subs(f,[xi],[x]));
fprintf('||∇f(xk)||: %d\t', norm(G)); 
fprintf('\n');

while norm(G) >= e
    %syms alpha; 
    
    x_old = x';
    for i=1:N 
      l(i)= x(i) + alpha * Pk(i);
    end
    
    for i=1:N
      x_f(i)= subs(f, xi, l); 
    end
    lh = subs(f,xi,x+alpha*Pk);
    rh = subs(f,xi,x)+c*alpha*G'*Pk;

    while lh>rh
        alpha = rho*alpha;
        lh = subs(f,xi,x+alpha*Pk);
        rh = subs(f,xi,x)+c*alpha*G'*Pk;
    end
    
    %GradF = diff(x_f, alpha);
    %alpha = solve(GradF, alpha); % Step size
    %alpha = 0.002
    for i=1:N 
      x(i)= x(i) + alpha * Pk(i); % Updated xk value
    end
   
    G = subs(df_dx, xi, x);
    Pk = -(G); % New Search Direction
    k = k+1;
    
    % display plot
    plot([x_old(1) x(1)],[x_old(2) x(2)],'*-r');
    % display table
    fprintf('k: %d\t', k); 
    for i=1:N
        fprintf('x%d: %d\t',i, x(i)); 
    end
    fprintf('f(xk): %d\t', double(subs(f,[xi], [x])));
    fprintf('||∇f(xk)||: %d\t', double(norm(G))); 
    fprintf('\n');
end
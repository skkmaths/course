clear all;
close all;
format long
%x = sym(zeros(1,2));
syms x [1,2]    % creates a vector of symbols [x1 x2]

x(1,1)   % prints the x1 symbol
x(1,2)   % prints the x2 symbol
x % print the symbol vector
f = x1^3 + x2^3;

A = subs(f,x,[0 1])     % computes f(0,1) 

% Gradient Computation:
for i=1:2 
    df_dx(i)= diff(f, "x"+i);
end

df_dx  % print the gradient vector

G =subs(df_dx, x, [0 1]) 

for i=1:2
    for j =1:2
        hf(i,j) = diff( diff (f, "x"+i), "x"+j);
    end
end

hf  % prints the hessian matrix 

hf11 = hf(1,1) % prints the first element at 1,1 of the hessian matrix

h = subs(hf,x,[0,1]) % computes the hessian at (0,1)




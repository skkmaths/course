% For arbitrary dimension use the symbol N

N = 3;
syms x [1,N];

f = x1^3 + x2^3+ x3^3;

for i=1:N 
    df_dx(i)= diff(f, "x"+i);
end

for i=1:N
    for j =1:N
        hf(i,j) = diff( diff (f, "x"+i), "x"+j);
    end
end

gradient = df_dx   % symbolic gradient vector
hessian = hf        % symbolic hessian matrix

G = subs(gradient, x, [0 1 1]) % eva

H = subs(hessian, x,[0 1 2])



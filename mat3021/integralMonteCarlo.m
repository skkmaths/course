a = 0; % lower limit of integration
b = 1; % upper limit of integration
n = 1000000; % number of random points to generate
sum = 0; % running total of function evaluations
f = @(s)  s.^2;
for i = 1:n
    x = a + (b-a) * rand(); % generate random point in [a,b]
    sum = sum + f(x); % add f(x) to running total
end
integral_estimate = (b-a) * sum / n; % compute estimate of integral
disp(integral_estimate);
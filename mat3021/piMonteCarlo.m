n = 10000000; % number of points to generate
count = 0; % number of points inside the circle
ind = zeros(1,n);
clearvars sum
for i = 1:n
    x = rand();
    y = rand();
    if x^2 + y^2 <= 1
        count = count + 1;
        ind(i) =1;
    end
    
end
pi_estimate1 = 4 * count / n;
pi_estimate2 = 4* sum(ind)/n;

disp(pi_estimate1);
disp(pi_estimate2);


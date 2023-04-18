n = 1000000; % number of points to generate
ind = zeros(1,n);
clearvars sum
for i = 1:n
    x = randi([1 1000]);
    if  mod(x,2)== 0 || mod(x,3) ==0 || mod(x,5)==0
        ind(i) =1;
    end
    
end

p = sum(ind)/n;

fprintf('Probability of getting a number divisible by 2, 3 or 5 is = %d.\n',p);
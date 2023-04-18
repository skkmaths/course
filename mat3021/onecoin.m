n = 1000000; % number of points to generate
ind = zeros(1,n);

clearvars sum
for i = 1:n
    x = randi([0 1]);
    if x == 1
        ind(i) =1;
    end
    
end

p = sum(ind)/n

fprintf('Probability of getting head %d \n= ', p);


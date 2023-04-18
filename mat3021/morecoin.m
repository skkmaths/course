n = 1000000; % number of points to generate
ind = zeros(1,n);
r =  4; %number of coin tossed
clearvars sum
for i = 1:n
    x = randi([0 1],1,r);
    if sum(x) == 3;
        ind(i) =1;
    end
    
end

p = sum(ind)/n;

fprintf('Probability of getting three heads = %d\n ', p);


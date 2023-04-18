N = 1000000; % number of points to generate
axis = [1:N];
clearvars sum
E = zeros(1,N);
S = 0;
for i = 1:N
    x = randi([1 6]);
    S = S+ x;
    E(i) = S/i;
end

fprintf('Expection = %d\n ', E(N));

%fplot(@x,[0 N])
hold on
plot(axis,E, 'LineWidth',2, axis)
xlabel('Trials');
ylabel('Expectation');
%xlim([0 1000]);
%ylim([-0.4 0.8]);

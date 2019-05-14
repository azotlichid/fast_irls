% TestM -- change in number of iterations in terms of m
rng(0);

M = [];

for i=1:30
    N = 200*i; % Signal length.
    T = 15;  % Signal sparsity.
    K = 150; % Number of measurements.
    
    x = zeros(N,1);
    q = randperm(N);
    x(q(1:T)) = sign(randn(T,1)); % Random +/- 1 signal.
    
    % Build measurement matrix.
    B = randn(K,N);
    B = orth(B');
    
    % Observations.
    d = B'*x;
    
    tic
    [r, f_energy, cnt1] = LinfSolve(B, d, .01, 0);
    t1 = toc;
        
    tic
    [r, f_energy, cnt2] = LinfSolve(B, d, .01, 1);
    t2 = toc;
    
    
    M = [M; N T K cnt1 t1 cnt2 t2];
end

save('test-m-linf.txt', 'M', '-ascii', '-tabs')

x = [1:size(M,1)];
h = plot(x, (M(:,4)), 'Color', 'red');
set(h,'LineWidth',3);
hold on;
h = plot(x, (M(:,6)),  'Color', 'blue');
set(h,'LineWidth',3);
axis tight;
set(h,'LineWidth',3);
saveas(h, 'test_m_linf_cnt.eps', 'epsc');
hold off

x = [1:size(M,1)];
h = plot(x, (M(:,5)), 'Color', 'red');
set(h,'LineWidth',3);
hold on;
h = plot(x, (M(:,7)),  'Color', 'blue');
set(h,'LineWidth',3);
axis tight;
set(h,'LineWidth',3);
saveas(h, 'test_m_linf_time.eps', 'epsc');
hold off


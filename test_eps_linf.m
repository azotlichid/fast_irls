% TestEps -- change in number of iterations in terms of eps
rng(0);

[B, d] = getprob(5);
M = [];

for i=1:12
eps = 1/2^i;
tic
[r, f_energy, cnt1] = LinfSolve(B, d, eps, 0);
t1 = toc    

tic
[r, f_energy, cnt2] = LinfSolve(B, d, eps, 1);
t2 = toc    

M = [M; eps cnt1 t1 cnt2 t2];
end

save('test-eps-linf.txt', 'M', '-ascii', '-tabs')

x = [1:size(M,1)];
h = plot(x, log(M(:,2)), 'Color', 'red');
set(h,'LineWidth',3);
hold on;
h = plot(x, log(M(:,4)),  'Color', 'blue');
set(h,'LineWidth',3);
axis tight;
set(h,'LineWidth',3);
saveas(h, 'test_eps_linf_cnt.eps', 'epsc');
hold off

x = [1:size(M,1)];
h = plot(x, (M(:,3)), 'Color', 'red');
set(h,'LineWidth',3);
hold on;
h = plot(x, (M(:,5)),  'Color', 'blue');
set(h,'LineWidth',3);
axis tight;
set(h,'LineWidth',3);
saveas(h, 'test_eps_linf_time.eps', 'epsc');
hold off
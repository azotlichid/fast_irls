rng(0);


[B,d] = getprob(5);

TEST = 0;

if TEST == 0
%tic
%fp = LinfSolveCvx(B, d);
%fprintf(' Cvx solver returned sol = %f\n', norm(fp, Inf));
%toc

tic
[r, f] = LinfSolve(B, d, 0.001, 1);
fprintf(' Solver returned sol = %f\n', norm(f, Inf));
toc
else
%tic
%f_cvx = LoneSolveCvx(B, d);
%fprintf('Cvx solver returned sol = %f\n', norm(f_cvx, 1));
%toc

fprintf('-----------\n');

tic
[r, f_energy] = LoneSolve(B, d, 0.01);
fprintf('Energy solver returned sol = %f\n', norm(f_energy, 1));
toc    

fprintf('-----------\n');

tic
f_candes = l1eq_pd(B*d, B', [], d, 0.01);
fprintf('Candes solver returned sol = %f\n', norm(f_candes, 1));
toc

fprintf('-----------\n');
end
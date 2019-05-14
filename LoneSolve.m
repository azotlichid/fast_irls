function [ c, f, cnt ] = LoneSolve( B, d, eps, use_line_search )
% Solve min ||f||_1 : B' * f = d
%   output: (1+eps) approximate solution consisting of
%           flow f and corresponding capacitances c
%           number of linear system solves cnt
if nargin >= 4
    LINE_SEARCH = use_line_search; % Set to 1 in order to attempt larger steps.
else
    LINE_SEARCH = 1;
end

[m, n] = size(B);
c = ones(m, 1); c = c / sum(c);
[x, E] = electric_pot(B, c, d);

cnt = 0;
best_ub = (d'*x)/norm(B*x, Inf);

while best_ub * (1+eps) < sqrt(E)
    cnt = cnt+1;
    OPT = sqrt(E)*(1-eps/2);
    alpha = next_alpha(B, c, d, OPT);
    
    if max(alpha) <= 1.000001
        fprintf('>> This is not supposed to happen << OPT = %f, sqrt(E) = %f\n', OPT, sqrt(E));
    end
    
    % Attempt to take longer step, in case it makes significant progress.
    if (sum(alpha-1) > 1) && LINE_SEARCH
        stride = 1;
        while 1
            cnt = cnt+1;
            cp = c.*(1+(2^stride*(alpha-1))); cp = cp / sum(cp);
            [~, Ep] = electric_pot(B, cp, d);
            if (Ep >= E) || (stride >= 20)
                stride = stride - 1;
                break;
            end
            stride = stride + 1;
        end
        c = c.*(1+(2^stride*(alpha-1))); c = c/sum(c);
    else
        c = c.* alpha;
        c = c / sum(c);
    end
    
    [x, E] = electric_pot(B, c, d);

    if (d'*x)/norm(B*x, Inf) > best_ub
        best_ub = (d'*x)/norm(B*x, Inf);
    end
end
fprintf('L1 solver finished after %d iterations\n', cnt);
f = diag(c)*B*x;
end

function [ alpha ] = next_alpha(B, c, d, OPT)
[x, E] = electric_pot(B, c, d);
alpha = ((B*x)/(d'*x)).^2  * OPT^2;
for i = 1:size(alpha)
    if alpha(i) < 1 /(1-eps)^2
        alpha(i) = 1;
    end
end
cc = c .* alpha;
end

function [ x, E ] = electric_pot(B, c, d)
Lap = B' * diag(c) * B;
Linvd = cgsolve(Lap, d, 1e-2, 10, 0); % Solve the linear system using conjugate gradient.
%Linvd = pinv(Lap)*d; % Solve the linear system using MATLAB's default methods.
x = Linvd;
E = d'*Linvd;
end
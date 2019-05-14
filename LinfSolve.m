function [ r, f, cnt ] = LinfSolve( B, d, eps, use_line_search )
%Solve min ||f||_inf : B' * f = d
%   output (1+eps) approximate solution consisting of
%          flow f and resistances r

if nargin >= 4
    LINE_SEARCH = use_line_search; % Set to 1 in order to attempt larger steps.
else
    LINE_SEARCH = 1;
end

[m, n] = size(B);
r = ones(m, 1); r = r / sum(r);
[f, E] = electric(B, r, d);

cnt = 0;

while E * (1+eps)^2 < (norm(f, inf))^2
    cnt = cnt + 1;
    OPT = sqrt(E) * (1+eps/2); 
    
    alpha = next_alpha(B, r, d, OPT);
    
    if (sum(alpha-1) > 1) && LINE_SEARCH
        stride = 1;
        while 1
            cnt = cnt+1;
            rp = r.*(1+(2^stride*(alpha-1))); rp = rp / sum(rp);
            [~, Ep] = electric(B, rp, d);
            if (Ep < E) || (stride >= 20)
                stride = stride - 1;
                break;
            end
            stride = stride + 1;
        end
        r = r.*(1+(2^stride*(alpha-1))); r = r/sum(r);
    else
        r = r.* alpha;
        r = r /sum(r);
    end
    [f, E] = electric(B, r, d);
end

fprintf('Linf solver finished after %d iterations\n', cnt);
end

function [ alpha ] = next_alpha(B, r, d, OPT)
[f, E] = electric(B, r, d);
m = size(B, 1);
alpha = (f/OPT).^2;
for i=1:size(alpha)
    if alpha(i) < 1+eps
        alpha(i) = 1;
    end
end
end

function [ f, E ] = electric(B, r, d)
Lap = B' * diag(1./r) * B;
Linvd =  cgsolve(Lap, d, 1e-2, 15, 0); %pinv(Lap)*d;
f = diag(1./r) * B * Linvd;
E = d'*Linvd;
end
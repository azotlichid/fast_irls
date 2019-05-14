function [ B, d ] = getprob(c)
if c == 0
    % Simple test example.
    B = [-1 0 0 1;
         -1 1 0 0;
         0 -1 1 0;
         0 0 -1 1;
         0 -1 0 1];
    d = [-1 0 0 1]';
elseif c == 1
    % Second simple test example.
    B = [-1  1  0 0;
          0 -1  1 0;
          0  0 -1 1;
         -1  0  1 0];
    d= [-1 0 0 1]';    
elseif c == 2
    % Random graph.
    B = random_graph(200, .1);
    [m, n] = size(B);
    d = rand(n,1); d = d - sum(d)/n;
elseif c == 3
    % Measurements based on sparse signal and random matrix.
    m = 900; n = 300;
    B = rand(m,n);
    x = zeros(m,1);
    for i = 1:round(m/20)
        x(randi(m)) = 1;
    end
    d = B'*x;
elseif c == 4
    N = 2000; % Signal length.
    T = 150;  % Signal sparsity.
    K = 1500; % Number of measurements.

    x = zeros(N,1);
    q = randperm(N);
    x(q(1:T)) = sign(randn(T,1)); % Random +/- 1 signal.

    % Build measurement matrix.
    B = randn(K,N);
    B = orth(B');
	
    % Observations.
    d = B'*x;
elseif c == 5
    N = 200; % Signal length.
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
else
    B = [0]; d = [0];
end
end

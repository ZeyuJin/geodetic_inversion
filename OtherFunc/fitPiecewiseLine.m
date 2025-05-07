function endpoints = fitPiecewiseLine(data, K)
%   endpoints = fitPiecewiseLine(data, K)
% Fitting fault segments coordinates given curved fault surface trace
% input
%   data :  N×2  matrix, first colun contitude, second column latitude
%   K    :  splines to fit (K ≥ 1 and K ≤ N‑1)
%
% output
%   endpoints : (K+1)×2  matrix，endpoint coordinates of each segment
%



validateattributes(data, {'numeric'}, {'2d','ncols',2});
x = data(:,1);  y = data(:,2);
N = numel(x);
K = min(max(round(K),1), N-1);


Sx  = [0; cumsum(x)];
Sy  = [0; cumsum(y)];
Sxx = [0; cumsum(x.^2)];
Sxy = [0; cumsum(x.*y)];
Syy = [0; cumsum(y.^2)];


err = inf(N,N);           
for i = 1:N
    for j = i:N
        n = j-i+1;
        sumx  = Sx(j+1)-Sx(i);
        sumy  = Sy(j+1)-Sy(i);
        sumxx = Sxx(j+1)-Sxx(i);
        sumxy = Sxy(j+1)-Sxy(i);
        sumyy = Syy(j+1)-Syy(i);

   
        denom = n*sumxx - sumx^2;
        if denom == 0   
            sse = sumyy - (sumy^2)/n;
        else
            a = (n*sumxy - sumx*sumy)/denom;      % slope
            b = (sumy - a*sumx)/n;                % intercept
            % SSE = Σ(y - (a x + b))²
            sse = sumyy + a^2*sumxx + n*b^2 - 2*(a*sumxy + b*sumy - a*b*sumx);
        end
        err(i,j) = sse;
    end
end


dp   = inf(K, N);  
bt   = zeros(K,N); 
dp(1,:) = err(1,:);

for k = 2:K
    for j = k:N
        best = inf; best_m = k-1;
        for m = k-1:j-1
            val = dp(k-1,m) + err(m+1,j);
            if val < best
                best = val; best_m = m;
            end
        end
        dp(k,j) = best;
        bt(k,j) = best_m;
    end
end


idx = zeros(K+1,1);
idx(end) = N;
for k = K:-1:1
    idx(k) = bt(k, idx(k+1));
end
idx(1) = 1;           


endpoints = zeros(K+1, 2);
for s = 1:K
    i = idx(s);   j = idx(s+1);
    x1 = x(i); x2 = x(j);
    n     = j-i+1;
    sumx  = Sx(j+1)-Sx(i);
    sumy  = Sy(j+1)-Sy(i);
    sumxx = Sxx(j+1)-Sxx(i);
    sumxy = Sxy(j+1)-Sxy(i);
    denom = n*sumxx - sumx^2;
    if denom == 0
        a = 0; b = sumy/n;    
    else
        a = (n*sumxy - sumx*sumy)/denom;
        b = (sumy - a*sumx)/n;
    end
    endpoints(s,  :) = [x1, a*x1 + b];
    endpoints(s+1,:) = [x2, a*x2 + b];
end
end
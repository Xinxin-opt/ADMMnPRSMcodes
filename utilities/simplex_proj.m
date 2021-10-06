function d = simplex_proj(c,tau)
%simplex_proj projects a vector c onto a simplex {x| x_i>=0, sum x_i = tau }

n = max(size(c));
p = -c;
pmax = max(p);
sm = sum(p);
if sm >= n*pmax - tau
    lambda = (tau+sm)/n;
    d = max(0, c + lambda);
    clear p;
    return;
end

p = sort(p);

sm = 0;
for i = 1:n-1
    smnew = sm + i*(p(i+1) - p(i));
    if smnew >= tau
        break
    end
    sm = smnew;
end

k = i;
delta = (tau - sm)/k;
lambda = p(k) + delta;
d = max(0, c + lambda);
%clear p;

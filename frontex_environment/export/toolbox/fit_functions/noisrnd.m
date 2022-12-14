function [x] = noisrnd(m,s)
% sample from noise distribution = truncated normal in [-.5,.5]
if numel(s) == 1
    s = s(ones(size(m)));
end
x = nan(size(m));
f = true(size(m));
while any(f(:))
    x(f) = normrnd(m(f),s(f));
    f(x >= -.5 & x <= .5) = false;
end
end
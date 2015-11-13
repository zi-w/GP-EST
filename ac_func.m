function val = ac_func(target, hyp, inf, mean, cov, lik, x, y, xs)
% AC_FUNC  the acquisition function used in the PI framework for GP
%   optimization. Intuitively, we would like to choose an x that gives us 
%   values close to the target value.

[yt,s2] = gp(hyp, inf, mean, cov, lik, x, y, xs);
val = (target - yt)./sqrt(s2);
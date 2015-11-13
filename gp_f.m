function [yy] = gp_f(xx,seed,prior)
% GP_F  A function sampled from a Gaussian process whose prior is 'prior'.
%   'seed' provides a random seed to this function.
%   GPML tool box is required.
%   For more information, see <a href="matlab: 
%   web('http://www.gaussianprocess.org/gpml/code/matlab/doc/index.html')">Documentation for GPML Matlab Code</a>


hyp = prior.hyp;
meanfunc = prior.meanfunc;
covfunc = prior.covfunc;
likfunc = prior.likfunc;

D = size(xx,2);

n = 100;

x = gpml_randn(seed, n, D);
K = feval(covfunc{:}, hyp.cov, x);
mu = feval(meanfunc{:}, hyp.mean, x);
y = chol(K)'*gpml_randn(1, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);

[yy] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, xx);

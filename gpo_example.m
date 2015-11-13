% This file is an example of Gaussian process (GP) optimization.
% It uses the EST algorithm, which is an extension of the PI frame work.
% The implementations of PI, and PI wiht known target value are also
% included in this framework.
%
% GPML tool box is required.
% For more information, see <a href="matlab: 
% web('http://www.gaussianprocess.org/gpml/code/matlab/doc/index.html')">Documentation for GPML Matlab Code</a>
%
% See also: gp_est.m


%% Specify the parameters of the problem
% dimension of input
D = 1;

% number of iterations
T = 100;

% the priors of the GP
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.1;1]; 
covfunc = {@covMaterniso,3}; ell = 0.1; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

prior.hyp = hyp;
prior.meanfunc = meanfunc;
prior.covfunc = covfunc;
prior.likfunc = likfunc;

% if isparafit = 1, the hyperparameters of GP will be optimized in each iteration. 
% otherwise it sticks with the priors specified above.
isparafit = 0;
prior.isparafit = isparafit;

% range of the search space for x
xlimit = [ones(1,D)*(-2);
    ones(1,D)*(2)];

% seed for the sampled function from the GP
seed = rand;

% the black box function we optimize. This can be replaced by any other
% function that the practitioner is optimizing. However, D and prior should
% be reset accordingly if this function is changed.
f = @(x)gp_f(x,seed,prior);

%% randomly generate the initial data
initx = rand_sample_x(1, D, xlimit);
inity = f(initx);

%%  prior.testid chooses which algorithm to use.
%   1 --- EST with numerical computation of m.
%   2 --- EST with approximate computation of m.
%   3 --- PI with epsilon shift (prior.pi) over current maximum.
%   4 --- PI with target value (prior.m).
prior.testid = 1;
% prior.pi = 0.1;

%% game start
% maxs stores all the target values in each iteration
maxs = zeros(1,T);
% xx stores all the inputs tested
xx = initx;
% yy stores all the output values from the inputs chosen
yy = inity;
for t = 1:T
    % choose the next point to evaluate
    [xnext, m] = gp_est(xx, yy, xlimit, prior);
    disp(['Round ' num2str(t) ': Estimated target value is m = ' num2str(m) '; Choose action x = ' num2str(xnext) ';']);
    % observe the output of the black box function (possibly with some noise)
    yt = f(xnext);% + normrnd(0,sn);
    disp(['f(' num2str(xnext) ') = ' num2str(yt)]);
    % update the current dataset
    xx = [xx;xnext];
    yy = [yy;yt];
    maxs(t) = m;
end


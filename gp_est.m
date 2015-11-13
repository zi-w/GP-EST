function [xnext, m] = gp_est(xx, yy, xlimit, prior)
% GP_EST  The EST algorithm for GP optimization under the framework of PI.
%   Theoretically EST works for any demension of x in a continuous
%   setting. This implementation uses a finite set of the space specified
%   by xlimit. xx and yy are the data collected, xlimit is the range of
%   search for the input, and prior is the prior of the GP we assumed.
%   The function outputs the next input to evaluate (xnext), and the
%   estimated maximum of the function (m).
%
%   prior.testid chooses which algorithm to use.
%   1 --- EST with numerical computation of m.
%   2 --- EST with approximate computation of m.
%   3 --- PI with epsilon shift (prior.pi) over current maximum.
%   4 --- PI with target value (prior.m).
%
%   GPML tool box is required.
%   For more information, see <a href="matlab:
%   web('http://www.gaussianprocess.org/gpml/code/matlab/doc/index.html')">Documentation for GPML Matlab Code</a>
%
%   See also: gpo_example.m

hyp = prior.hyp;
meanfunc = prior.meanfunc;
covfunc = prior.covfunc;
likfunc = prior.likfunc;
testid = prior.testid;
isparafit = prior.isparafit;

% optimize over the hyperparameters
if isparafit
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, xx, yy);
end

% dimension of the input
D = size(xx,2);

m = estm;

initx = rand_sample_x(1, D, xlimit);
% use unconstrained pattern search minimization to minimize ac_func
af = @(x)ac_func(m,hyp, @infExact, meanfunc, covfunc, likfunc, xx, yy, x);
xnext = patternsearch(af, initx,[],[],[],[],xlimit(1,:),xlimit(2,:));

    function m = estm
        % estimate the max of the function
        if testid == 1 % numerical integratation
            % random sample RS points and evaluate the mean and variance
            RS = 10000;
            xt = rand_sample_x(RS, D, xlimit);
            [yt,s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xx, yy, xt);
            s2 = sqrt(s2);
            % calculate m
            m = max(yy);
            m0 = m;
            intcnt = 1;
            logprodphi = zeros(100,1);
            while intcnt == 1 || logprodphi(intcnt-1) < 0
                logprodphi(intcnt) = sum(log(normcdf((m0-yt)./s2, 0, 1)))';
                m = m + (1-exp(logprodphi(intcnt)))*0.05;
                m0 = m0+0.05;
                intcnt = intcnt + 1;
            end
            
        elseif testid == 2 % approximate integration with Gaussian
            % random sample RS points and evaluate the mean and variance
            RS = 10000;
            xt = rand_sample_x(RS, D, xlimit);
            [yt,s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xx, yy, xt);
            s2 = sqrt(s2);
            % calculate m
            m = max(yy);
            m0 = m;
            logprodphi10 = sum(log(normcdf((m0+1-yt)./s2, 0, 1)));
            if log(1-exp(logprodphi10)) == 0
                tempsig = 50;
            else
                tempsig = sqrt(-1/2/log(1-exp(logprodphi10)));
            end
            logprodphi1 = sum(log(normcdf((m0-yt)./s2, 0, 1)));
            m = m + 1/normpdf(0,0,tempsig)*(1-exp(logprodphi1))/2;
        elseif testid == 3 % pi
            if isfield(prior,'pi')
                m = max(yy)+prior.pi;
            else
                disp('Error: prior.pi needs to be specified.');
            end
        elseif testid == 4 %
            if isfield(prior,'m')
                m = prior.m;
            else
                disp('Error: prior.m needs to be specified.');
            end
        else
            disp('Error: unknown testid.');
        end
    end

end
function fit = fit_baseline(theta, count, useBootstrapping)
% fit = fit_baseline(theta, R)
% fit baseline with proper min and max using maximum likelihood with
% poisson distributed spike counts
% Assumes directions are in degrees
if nargin < 3
    useBootstrapping = false;
end

% initialize parameters
minFR = min(count); 
maxFR = max(count); 
params0 = [0.5*(minFR+maxFR)];
if (params0 == 0)
    params0 = 0.001;
    maxFR = 1;
end
opts = optimset('MaxFunEval', 10e3, 'Display', 'off');

% if you want to try bounds
LB = [0.00001]; % lower bound
UB = [maxFR];

% LB = [];
% UB = [];
if useBootstrapping
    nBoot = 100; %2000;
    pBoot = zeros(nBoot,numel(params0));
    
    nTrials = numel(count);
    
    for i = 1:nBoot
        inds = randi(nTrials, nTrials, 1);
        fun = @(params) VonMises.neglogli_poissGLM(VonMises.baseline(theta,params), count(inds));
        pBoot(i,:) = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
    end
    
    paramsSD = prctile(pBoot, [16 84]); % 68% confidence intervals
    
    %build objective function
    fun = @(params) VonMises.neglogli_poissGLM(VonMises.baseline(theta,params), count);
    
    % optimization
    [phat, fval, ~, ~, ~, ~, H] = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
else
    %build objective function
    fun = @(params) VonMises.neglogli_poissGLM(VonMises.baseline(theta,params), count);
    
    % optimization
    [phat, fval, ~, ~, ~, ~, H] = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
    
    % error bars
    paramsSD = sqrt(diag(inv(H)))';
end

% fit output
fit.paramsML = phat;
fit.fvalue = fval;
fit.paramsSD = paramsSD;
if exist('pBoot', 'var')
    fit.pBoot = pBoot;
    fit.FR = phat(1);
    fit.FRSD = abs(prctile(pBoot(:,1), cibnd) - phat(1));
else
    fit.FR = phat(1);
    fit.FRSD = paramsSD(1);
end


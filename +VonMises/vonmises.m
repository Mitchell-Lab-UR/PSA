function y = vonmises(x, params)
%VONMISES evaluate von mises function, inputs in RADIANS
%
% y = VONMISES(x, params) evaluates at x using the given params, where
% params(1) is min val, params(2) is max val, params(3) is kappa (width
% parameter) and params(4) is preferred orientation
%
% NOTE that this formulation of the von mises function removes some
% 'couplings' between terms that would otherwise make optimizing it
% difficult. For example, min and max are used rather than min and
% amplitude. This way the max val is easy to bound at 1.0 for bernoulli
% likelihoods. Notice also that the exp(cos()) term has likewise been normalized
% to range in [0,1]. This is easy to derive noting that cos() is in [-1,1]
% and plugging in those bounds to 'remap' it
%
% ranges [exp(-k), exp(k)]
% now subtract min, divide by max-min to get to [0,1]
%
% (exp(k*cos(...)) - exp(-k)) / (exp(k) - exp(-k))
% now multiply by exp(+k) in all 4 terms
%
% (exp(k*cos(...))exp(k) - 1) / (exp(2k) - 1)
% (exp(k*(cos(...)+1)) - 1) / (exp(2*k) - 1)

base = params(1);
amp = params(2);
%minval = params(1);
%maxval = params(2);
k = params(3);
pref = params(4);

%****** I don't follow the issue for base and amp parameters vs min and max
%****** the bound at max is already limited by maxFR*1.5
if (k > 0)
   y = base  + amp * (exp(k*(cos(x - pref)+1)) - 1) / (exp(2*k) - 1);  
else   % modified here below to allow full range on kappa, broader peaks
    if (k == 0)
       y = base + amp * (0.5 * (cos(x - pref) + 1));  
    else   % for k negative, flip function to make broad at the peak
       y = base  + amp * (1-(exp(-k*(cos(x - pref + pi)+1)) - 1) / (exp(-2*k) - 1) );      
    end
end

% y = minval + (maxval - minval) * (exp(k*(cos(x - pref)+1)) - 1) / (exp(2*k) - 1);
% y = minval + (maxval - minval) * (exp(k*(cos(x - pref)+1)) - 1) / (exp(2*k) - 1);
end
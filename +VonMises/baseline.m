function y = baseline(x, params)
%BASELINE evaluate baseline function, inputs in RADIANS
%
% y = BASELINE(x, params) evaluates at x using the given params, 
%     and simply fitting same baseline for all values
%

y = ones(size(x))*params(1);

end
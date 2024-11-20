function info = Compute_mutualInformation(x, y)
% calculates the mutual information (bits) between two vectors of data.
%
%
% Author - Hayden Scott
% June 2020

% I(X;Y) = D_kl(P(x,y) || Px xo Py ) 
% mutual information is related to the kl divergence
%% Sanity checks:
assert(numel(x)==numel(y), 'X and Y must have paired observations');
% force column vectors
x = x(:);
y = y(:);

% if only one observed unique data point 
if numel(unique(y))==1||numel(unique(x))==1
    info = 0;
    return
end

%% Bins
% this function provides an unbiased binning of the data

function bins = goodBins(dat)
    % creates unbiased bins using the formula from Scott 1979
    % this is the same formula built-in matlab functions use
    binWidth = 3.49*std(dat)*numel(dat)^(-1/3);
    bins = min(dat)-binWidth:binWidth:max(dat)+binWidth;
end

xBins = goodBins(x);
yBins = goodBins(y);

%% Calculate marginal distributions
% probability of observing the data in x and y independent of each other
[margx, ~, inBinX] = histcounts(x,xBins);
margx = margx./sum(margx);

[margy, ~, inBinY] = histcounts(y,yBins);
margy = margy./sum(margy);

%% Calculate the joint distribution
% probability of observing data in y given x
% the same as observing data in x given y

joint = accumarray([inBinX,inBinY],1,[numel(xBins)-1,numel(yBins)-1],@sum);
joint = joint./sum(joint(:));

%% calculate Mutual Information

marges = repmat(margx,1,size(margy,2)).*repelem(margy,1,size(margx,2));
val = joint(:).*log2(joint(:)./(marges'));
info = nansum(val(val~=Inf));

end

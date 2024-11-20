function smo = gauss_smooth(psth,Gsig)

    % Make the number of samples depending on the gaussian window size
    gaussian_filter_size = 4*Gsig-1; % if Gsig = 10, 19 samples total
                                     % 9 left & 9 right from the mean

    % Make smoothing kernel using gaussian filter
    for i = 1:gaussian_filter_size
        gauss  = exp(-(((i-(2*Gsig)).^2)./(2*Gsig^2)));
        gauss_filter(i,:) = gauss;
    end
    % Normalize the gaussian filter
    gauss_smooth = gauss_filter/sum(gauss_filter);
    psth_size    = length(psth);
    filter_size  = length(gauss_smooth);
    filter_cent = floor((filter_size+1)/2);

    for i=1:psth_size   % size_smooth

        % Always 0 for the initial value (only sum from product of two vectors)
        smo(i) = 0;
        nomo(i) = 0;

        % Apply filter to data
        for j = 1:filter_size
             diff = (j-filter_cent);   % this coordinate, goes - 3*sig up to 3*sig
             samp = (i+diff);
             if ( (samp >= 1) && (samp <= psth_size) )
                 if (j == i) && isnan(psth(samp))  % if center val Nan, keep Nan
                     smo(i) = NaN;
                 else
                   if ~isnan(psth(samp))
                     smo(i) = smo(i) + (psth(samp) * gauss_smooth(j));
                     nomo(i) = nomo(i) + gauss_smooth(j);
                   end
                 end
             end        
        end
        %********
        if (nomo(i) > 0)
            smo(i) = smo(i) / nomo(i);
        end
        %***********
    end
    
return;
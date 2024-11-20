function [tt,uu,su] = CompPSTH(StimRast,StimWin,N,Smoothing)
% function [tt,uu,su] = CompPSTH(StimRast,StimWin,N,Smoothing)
%*** inputs:   StimRast - Nx2 where col 1 is time and col 2 is trials
%**            StimWin - Time interval to put into raster
%**            N - Number of total trials put into raster
%**            Smoothing - use Gaussian smoothing of this sigma
%**
%*** outputs:
%**         tt = time vec in secs
%**         u = mean PSTH
%**         su = sem of PSTH

   %***** build matrix form of raster with 1ms bins
   T = floor(1000*(StimWin(2)-StimWin(1))) + 1;
   Raster = zeros(N,T);
   if ~isempty(StimRast)
     ti = ceil( (StimRast(:,1)-StimWin(1))*1000);
     ni = StimRast(:,2);
     for k = 1:length(ti)
       if (ni(k) > 0) && (ti(k) > 0)
          Raster(ni(k),ti(k)) = 1;
       end
     end
   end
   
   %******* compute the Jacknife
   smoothsub = [];
   JNum = 10;   % the number of Jacknife's in estimate
   for i = 1:JNum 
     aex = 1+floor((i-1)*N/JNum);
     bex = ceil(i*N/JNum);
     excludeset = [1:aex,bex:N];
     traces = mean(Raster(excludeset,:))*1000;
     smooth_data2 = spikestats.gauss_smooth(traces,Smoothing);
     smoothsub = [smoothsub; smooth_data2];
   end
   %***********
   tt = ((1:T)*(1/1000)) + StimWin(1);
   uu = mean(smoothsub);
   su = std(smoothsub) * sqrt(JNum-1);
      
return;

% 
% function smo = gauss_smooth(psth,Gsig)
% 
%     % Make the number of samples depending on the gaussian window size
%     gaussian_filter_size = 4*Gsig-1; % if Gsig = 10, 19 samples total
%                                      % 9 left & 9 right from the mean
% 
%     % Make smoothing kernel using gaussian filter
%     for i = 1:gaussian_filter_size
%         gauss  = exp(-(((i-(2*Gsig)).^2)./(2*Gsig^2)));
%         gauss_filter(i,:) = gauss;
%     end
%     % Normalize the gaussian filter
%     gauss_smooth = gauss_filter/sum(gauss_filter);
%     psth_size    = length(psth);
%     filter_size  = length(gauss_smooth);
%     filter_cent = floor((filter_size+1)/2);
% 
%     for i=1:psth_size   % size_smooth
% 
%         % Always 0 for the initial value (only sum from product of two vectors)
%         smo(i) = 0;
%         nomo(i) = 0;
% 
%         % Apply filter to data
%         for j = 1:filter_size
%              diff = (j-filter_cent);   % this coordinate, goes - 3*sig up to 3*sig
%              samp = (i+diff);
%              if ( (samp >= 1) && (samp <= psth_size) )
%                  smo(i) = smo(i) + (psth(samp) * gauss_smooth(j));
%                  nomo(i) = nomo(i) + gauss_smooth(j);
%              end       
%         end
%         %********
%         if (nomo(i) > 0)
%             smo(i) = smo(i) / nomo(i);
%         end
%         %***********
%     end
%     
% return;
% 

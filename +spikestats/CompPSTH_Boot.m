function [tt,uu,su] = CompPSTH_Boot(StimRast,StimWin,N,M,Smoothing)
% function [tt,uu,su] = CompPSTH_Boot(StimRast,StimWin,N,M,Smoothing)
%*** inputs:   StimRast - Nx2 where col 1 is time and col 2 is trials
%**            StimWin - Time interval to put into raster
%**            N - Number of total trials put into raster
%**            M - Number of replicated distributions
%**            Smoothing - use Gaussian smoothing of this sigma
%**
%*** outputs:
%**         tt = time vec in secs
%**         u = mean PSTH
%**         su = sem of PSTH

   %**** one issue is I've over-duplicated sets, which shrinks the
   %**** error bars.
   % *** INSTEAD, get error bars for each set, then average them
   %*********

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
   
   NN = (N/M); % number per replicate
   tts = [];
   uus = [];
   sus = [];
   for kk = 1:M   
     IRaster = Raster((1+(kk-1)*NN):(kk*NN),:);  
     %******* compute the Jacknife
     smoothsub = [];
     JNum = 10;   % the number of Jacknife's in estimate
     for i = 1:JNum 
       aex = 1+floor((i-1)*NN/JNum);
       bex = ceil(i*NN/JNum);
       excludeset = [1:aex,bex:NN];
       traces = mean(IRaster(excludeset,:))*1000;
       smooth_data2 = spikestats.gauss_smooth(traces,Smoothing);
       smoothsub = [smoothsub; smooth_data2];
     end
     %***********
     tt = ((1:T)*(1/1000)) + StimWin(1);
     uu = mean(smoothsub);
     su = std(smoothsub) * sqrt(JNum-1);
     %********
     tts = [tts ; tt];
     uus = [uus ; uu];
     sus = [sus ; su];
     %***
   end
   tt = mean(tts);
   uu = mean(uus);
   su = mean(sus);
   
return;


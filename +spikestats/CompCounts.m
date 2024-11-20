function OriSpk = CompCounts(ta,tb,OriInd,StimRast)
%  function OriSpk = CompCounts(ta,tb,StimRast)
%   inputs:  ta - start time interval
%            tb - end time interval
%            OriInd - indexes of orientation values per trial
%            StimRast - raster of spike times (list of trial by spike time)
%                  first column is times, second is trials

       %***** get counts per trial in window
       OriSpk = zeros(size(OriInd));
       if (isempty(StimRast))
           return;
       end
       for k = 1:length(OriInd)
           kz = find( StimRast(:,2) == k);  % find trial
           sp = find( (StimRast(kz,1) >= ta) & (StimRast(kz,1) < tb) );
           if ~isempty(sp)
               rate = length(sp)/(tb-ta);
           else
               rate = 0;
           end     
           OriSpk(k) = rate;
       end

 return;
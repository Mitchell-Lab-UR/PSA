function [tt,uu,su] = CompCountSlider(StimRast,PrefOri,NPrefOri,OriVals,OriInd,LockWin,CountWin,Type)
%  function [tt,uu,su] = CompCountSlider(StimRast,PrefOri,NPrefOri,OriVals,OriInd,LockWin,CountWin)
%*** inputs:   StimRast - Nx2 where col 1 is time and col 2 is trial ind
%**            POri - preferred ori (in integer)
%**            OriVals - list of orientations, in degs
%**            OriInd - integer of orientation per trial list
%**            LockWin - total time window to compute over with sliding win
%**            CountWin - counting window
%**            Type - (1, AUC), (2, Fano), (3, R Ori)
%*** outputs:
%**         tt = time vec in secs
%**         u = mean AUC
%**         su = sem of AUC

   %***** loop over time with the counting window
   Step = 0.005;  % 5 ms steps
   tt = [];
   uu = [];
   su = [];
   OriSpk = zeros(size(OriInd));
   for t = LockWin(1):Step:LockWin(2)
       tt = [tt t];
       %***********
       ta = max(LockWin(1),(t-(CountWin/2)));
       tb = min(LockWin(2),(t+(CountWin/2)));
       OriSpk = spikestats.CompCounts(ta,tb,OriInd,StimRast);  % get counts from Rast in the window [ta,tb]
       %***** get counts per trial in window
       if (Type == 1)
           [uval,sval,ci] = spikestats.Compute_AUC(PrefOri,NPrefOri,OriInd,OriSpk);
       end
       if (Type == 2)
           
       end
       if (Type == 3)
           [uval,sval] = Compute_Resultant(OriVals,OriInd,OriSpk);
       end
       if (Type == 5)
           [uval,sval,ci] = spikestats.Compute_AUC_Circle(OriVals,OriInd,OriSpk);
       end
       uu = [uu uval];
       su = [su sval];
       %******
   end   
return;


function [uval,sval] = Compute_Resultant(OriVals,OriInd,OriSpk);
       
        %*** Jacknife and compute AUC
        smoothsub = [];
        JNum = 10;   % the number of Jacknife's in estimate
        N = length(OriInd); 
        for i = 1:JNum 
               aex = 1+floor((i-1)*N/JNum);
               bex = ceil(i*N/JNum);
               eset = [1:aex,bex:N];
               %**********
               cval = 0 + j * 0 ;
               wval = 0;
               for k = 1:length(eset)
                   kk = eset(k);
                   cval = cval + OriSpk(kk) * exp( j * (OriVals( OriInd(kk)) * (pi/180)) );
                   wval = wval + OriSpk(kk);
               end
               if (wval > 0)
                   val = abs(cval) / wval;
               else
                   val = 0;
               end
               %**********
               smoothsub = [smoothsub; val];
        end
        uval = nanmean(smoothsub);
        sval = nanstd(smoothsub) * sqrt(JNum-1);
        %******
return;



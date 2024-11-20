function [val,sval,ci] = Compute_AUC_Circle(OriVals,OriInd,OriSpk)
%   function [val,sval] = Compute_AUC_Circle(OriInd,OriSpk)
%   input:  OriVals - list of possible values
%           OriInd - index of orientation integers per trial
%           OriSpk - index of spike rates per trial
%
%  computes the AUC between every 3 adjacent indexes      
%   ouput:  val - AUC value, sem, and conf intervals (ci two element)

      N = length(OriVals);
      %*******
      uu = [];
      su = [];
      ci = [];
      %*******
      
      for kp = 1:N
         
         %********* 
         if (kp == 1)
             PrefOri = [N,1,2];
         else
             if (kp == N)
                 PrefOri = [(N-1),N,1];
             else
                 PrefOri = [(kp-1),kp,(kp+1)];
             end
         end
         %******
         NPrefOri = PrefOri - floor(N/2);
         zz = find(NPrefOri < 1);
         NPrefOri(zz) = NPrefOri(zz) + N;
         %**********
         
         %******* list of pref counts
         alist = [];
         for k = 1:length(PrefOri)
            z = find( OriInd == PrefOri(k));
            if ~isempty(z)
              alist = [alist ; OriSpk(z)];
            end
         end
         %****** list of non pref counts
         blist = [];
         for k = 1:length(NPrefOri)
            z = find( OriInd == NPrefOri(k));
            if ~isempty(z)
              blist = [blist ; OriSpk(z)];
            end
         end
         %**** must define alist as high rate *******
         if mean(blist) > mean(alist)
             tmp = alist;
             alist = blist;
             blist = tmp;
         end
         %************
         if isempty(alist) | isempty(blist)
            val = NaN;
            sval = NaN;
            ci = [NaN,NaN];
         else
            auc = spikestats.JakeROC([alist ; blist]',...
                             [ones(size(alist)) ; zeros(size(blist))]',0.05);
            val = auc.AUC;
            sval = 0.25*(max(auc.ci)-min(auc.ci));
            ci = auc.ci;
         end
         %*************
         uu = [uu ; val];
         su = [su ; sval];
         ci = [ci ; ci];
      end
      %******** average the results across samples
      val = nanmean(uu);
      sval = nanmean(su)/sqrt(N/3);  % number of indep comparisons
      ci = nanmean(ci);
      %**********
      
return


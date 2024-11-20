function [val,sval,ci] = Compute_AUC(PrefOri,NPrefOri,OriInd,OriSpk)
%   function [val,sval] = Compute_AUC(PrefOri,NPrefOri,OriInd,OriSpk)
%   input:  PrefOri  - list of integers for pref ori (related to OriInd
%           NPrefOri - list of integers for non-pref ori
%           OriInd - index of orientation integers per trial
%           OriSpk - index of spike rates per trial
%
%  computes the AUC between PrefOri and NprefOri distributions     
%   ouput:  val - AUC value, sem, and conf intervals (ci two element)

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
      
return


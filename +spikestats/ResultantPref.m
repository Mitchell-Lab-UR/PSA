function kpref = ResultantPref(OriVals,Wvec)
%
% function Pref = ResultantPref(OriVals,Wvec)
% inputs:  OriVals:  1xN vector of orientation values (in degs)
%          Wvec:  vector to align with best Orival and return integer
      
      kpref = 0;
      dotprod = -Inf;
      for ko = 1:length(OriVals)
          vc = exp( i * ((OriVals(ko)*(pi/180))));
          dp = sum( (real(Wvec)*real(vc)) + (imag(Wvec)*imag(vc)) );
          if (dp > dotprod)
            kpref = ko;
            dotprod = dp;
          end
      end
      
return


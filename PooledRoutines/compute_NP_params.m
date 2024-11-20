% non-parametric estimates of tuning parameters
function [ampa,ampu,basea,baseu,awid,uwid] = compute_NP_params(ratt,rutt)
          %**** assumes rank order around mean preference already done
%           basea = nanmean(ratt(1:4));
%           baseu = nanmean(rutt(1:4));
%           ampa = nanmean(ratt(13:16))-basea;
%           ampu = nanmean(rutt(13:16))-baseu;  
%           zaa = max((ratt-basea),0);
%           zaa = min((zaa/ampa),1);  % bound from 0 to 1
%           zuu = max((rutt-baseu),0);
%           zuu = min((zuu/ampu),1);  % bound from 0 to 1
%           awid = sum(zaa)/length(zaa);
%           uwid = sum(zuu)/length(zuu);
        
          basea = nanmean(ratt(1:4));
          baseu = nanmean(rutt(1:4));
          ampa = nanmean(ratt(13:16))-basea;
          ampu = nanmean(rutt(13:16))-baseu;  
          zaa = max((ratt-basea),0);
          zaa = min((zaa/ampa),1);  % bound from 0 to 1
          zuu = max((rutt-baseu),0);
          zuu = min((zuu/ampu),1);  % bound from 0 to 1
          %******* find mid-point in mass
          tot = sum(zaa);  % total
          otso = 0;
          tso = 0;
          mid = NaN;
          for i = 1:length(zaa)
              tso = (sum(zaa(1:i))/tot);
              if (tso > 0.5)
                 mid = ( i*(0.5-otso) + (i-1)*(tso-0.5) )/(tso-otso);
                 break;
              else
                 otso = tso; 
              end
          end
          awid = length(zaa)-mid;
          %********
          tot = sum(zuu);  % total
          otso = 0;
          tso = 0;
          mid = NaN;
          for i = 1:length(zuu)
              tso = sum(zuu(1:i))/tot;
              if (tso > 0.5)
                 mid = ( i*(0.5-otso) + (i-1)*(tso-0.5) )/(tso-otso);
                 break;
              else
                 otso = tso;
              end
          end
          uwid = length(zuu)-mid;   % from peak down to 0.5
          %*******  
return;
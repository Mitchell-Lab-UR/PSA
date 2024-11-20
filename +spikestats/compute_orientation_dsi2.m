%****** compute the orientation tuning and DSI value
function OriTune = compute_orientation_dsi2(OriVals,OriInd,OriSpk)
  % function OriTune = compute_orientation_dsi(OriVals,OriInd,OriSpk)
  %  *** Inputs:  OriVals:  list of Orientation values (in degs)
  %               OriInd:   integer index per trial into OriVals
  %               OriSpk:   spike rate from that trial
  %  *** Output:  Various circular statistics and tuning curve with SEM
  
  %******** compute the mean and sem of rate per orientation
  NVals = length(OriVals);
  OriTune.OriVals = OriVals;
  OriTune.Mu = NaN(1,NVals);
  OriTune.Num = NaN(1,NVals);
  OriTune.Sem = NaN(1,NVals);
  for k = 1:NVals
       z = find( OriInd == OriVals(k));
       if ~isempty(z)
          OriTune.Mu(k) = nanmean(  OriSpk(z) );
          OriTune.Num(k) = length(z);
          if (OriTune.Num(k))
             OriTune.Sem(k) = nanstd( OriSpk(z) ) / sqrt(OriTune.Num(k));
          end
       end
  end
  %OriVals
  %OriTune.Mu
  %OriTune.Num
  %OriTune.Sem
  %input('check');
  %********* then compute the DSI and preferred vector
  % but note, is done from raw mean values, not Von Mises fit
  wsum = 0;
  Wvec = 0;
  for k = 1:NVals
     ori = OriVals(k)*(pi/180);
     vc = complex(cos(ori),sin(ori));
     if ~isnan(OriTune.Mu(k))
        wsum = wsum + OriTune.Mu(k); 
        Wvec = Wvec + (OriTune.Mu(k) * vc);
     end
  end
  if (wsum > 0)
      OriTune.DSI = abs( Wvec / wsum );  % note DSI is r in Berens circ_var
      OriTune.AngDev = sqrt( 2*(1-OriTune.DSI));   
      OriTune.AngStd = sqrt( -2 * log(OriTune.DSI) );
      OriTune.Wvec = Wvec / wsum;
      %***** search to find the preferred (integer) orientation
      OriTune.Pref = spikestats.ResultantPref(OriVals,OriTune.Wvec);
  else
      OriTune.DSI = 0;
      OriTune.AngDev = sqrt(2);    
      OriTune.AngStd = Inf;
      OriTune.Wvec = [0 + 0*j];
      OriTune.Pref = NaN;  % no preferred value
  end
  
  return;


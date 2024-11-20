function [Afit,Atune,Awid,Araw,Fval,BFval,BFrate] = FitWithVonMises(OriVals,JackN,Ori,Spk,OFit)
%***** function [Afit,Atune,Awid,Araw,Fval] = FitWithVonMises(OriVals,JackN,Ori,Spk)
%*** inputs:  OriVals - direction vals(0-360) of discrete sampled directions (usually 16 of them)
%***          JackN - number Jacknifes to compute confidence intervals
%***          Ori - Nx1 vector of orientation values, where N is number of trials
%***          Spk - Nx1 vector of spike counts on each trial
%***          OFit - 1x4 param vector, start point in fits
%****************

  %*** possible that might call routine with no spikes, or no trials
  if isempty(Spk) || (sum(Spk)==0)
      Afit = [];  % look for empty return to throw out unit
      Atune = [];
      Awid = [];
      Araw = [];
      Fval = NaN;
      BFval = NaN;
      BFrate = NaN;
      return;
  end
  
  %***********
  NA = size(Spk,1);
  mfitA = [];
  tuneA = [];
  wfitA = [];
  fvalA = [];
  %****** fit a baseline 
  fitB = VonMises.fit_baseline(Ori,Spk,0);
  BFval = fitB.fvalue;
  BFrate = fitB.FR;
  %******
  MissedJack = 0;
  for ik = 1:JackN
      ia = 1 + (ik-1)*floor(NA/JackN);
      ib = ik * floor(NA/JackN);
      subset = [1:ia, ib:NA];
      %****** Call to VonMises function fitting
      fitA = VonMises.fit_vonmises(Ori(subset), Spk(subset), 0, OFit);  % not using Bootstrap ... here Jacknife instead
      if isempty(fitA)
         MissedJack = MissedJack + 1;   % selected a subset that had no spikes (must be low rate unit!)
         continue;
      end
      %*********
      for k = 1:size(OriVals,1)
         tuna(k) = VonMises.vonmises(OriVals(k)*(pi/180),fitA.paramsML);
      end
      tuneA = [tuneA ; tuna];
      mfitA = [mfitA ; fitA.paramsML];
      fvalA = [fvalA ; fitA.fvalue];
      %******* compute half-width
      for k = 1:360
          ztuna(k) = VonMises.vonmises((k*(pi/180)),fitA.paramsML);
      end
      maxo = max(ztuna);
      mino = min(ztuna);
      half = 0.5*(maxo+mino);
      zz = find( ztuna > half);
      wid = length(zz);
      %**************************
      wfitA = [wfitA ; wid];
  end
  %******
  if (MissedJack)
      JackN = JackN - MissedJack;
      if (JackN < 2)
          disp('Failed to fit Von Mises Entirely .... low rate??');
          Afit = [];
          Awid = [];
          Araw = [];
          return;
      end
  end
  %********
  mtuneA = mean(tuneA);
  stuneA = std(tuneA) * sqrt(JackN-1);
  Atune.mu = mtuneA;
  Atune.sem = stuneA;
  %******
  mfitmu = mean(mfitA);
  mfitsem = std(mfitA) * sqrt(JackN-1);
  Afit.mu = mfitmu; 
  Afit.sem = mfitsem; 
  %***** last parameter is an angle, so we need to convert coordinates
  ango = 0;
  for i = 1:JackN
      ango = ango + (mfitA(i,2) * exp( j * mfitA(i,4)));
  end
  amu = angle(ango);
  if (amu < 0)
      amu = amu + (2*pi);
  end
  Afit.mu(4) = amu;
  varo = [];
  for i = 1:JackN
       dang = mfitA(i,4) - amu;
       if (dang > pi)
         dang = dang - (2*pi);
       end
       if (dang < -pi)
         dang = dang + (2*pi);
       end
       varo = [varo ; dang];
  end
  Afit.sem(4) = std(varo) * sqrt( JackN - 1);
  %*********
  Fval = mean(fvalA);  % return goodness of fit
  %*******
  Awid.mu = mean(wfitA);
  Awid.sem = std(wfitA) * sqrt( JackN - 1);
  %******* compute the raw tuning curve as well
  Araw = [];
  for k = 1:size(OriVals,1)
      zz = find( Ori == OriVals(k) );
      Araw(k) = mean( Spk(zz));
  end
  %************************
  
 return;
  

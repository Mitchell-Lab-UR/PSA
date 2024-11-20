function Info = Compute_DSI(Info,Exp,Unit,H)
%******* 
%******  function Info = Compute_DSI(Info,Exp,Unit,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          Unit - spike cluster to run computation on
%***          H - handle to plot window if provided, or if 1 create, else []
%***
%*** Outputs: Info - it computes the DSI from all trial types,
%***                   and puts confidence intervals on it

  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh = figure('position',[100 100 1200 400]); 
       HA = subplot('position',[0.0752 0.15 0.25 0.7]);
       HB = subplot('position',[0.4 0.15 0.25 0.7]);
       HC = subplot('position',[0.725 0.15 0.25 0.7]);
    else
        if isempty(H)
          return;
        else
          if (length(H) ~= 2)
              disp('Error with subplot pass to Compute_DSI, requires two subplots');
              return;
          else
              HA = H{1};
              HB = H{2};
              HC = H{3};
          end
        end
    end
    %*******
    plot_compute_dsi(Info,HA);
    plot_tuning(Info,HB);
    plot_stimpsth(Info,HC);
    %********
    if isfloat(H)  && (H == 2)
       z = getframe(hh);  % current figure
       uname = [Info.pathplot,filesep,'DSI_',Info.tagname,'.png'];
       disp(sprintf('Storing image of DSI graph at %s',uname)); 
       imwrite(z.cdata,uname); % store image for later review
       close(hh);
    end   
    return;
  end
  
  %****** parameters for computing DSI and neuron inclusion
  DSI_threshold = 0.05;   % flag to include unit or not, based on the DSI being
                          % significantly above this threshold with p=0.05
  JackN = 10;             % Jacknife to get confidence intervals
  StimInt = [0.04,0.20];  % Interval to compute motion selectivity (DSI), % [0.04,0.14]
                          % note, used this wider interval in original analysis
  StimWin = [-0.1,0.25];   % plot stim locked firing rate PSTH 
  Smooth = 5;             % Gaussian sig for smooth PSTH plots
  %******
  THETAERROR = 90;        % Require Von Mises Pref Theta SEM less than this (not strict)
                          % We can filter unit based on R2 instead later
                          % from the Von Mises fit, so no need to set limit
                          % here
  %*******************
  
  %****** sanity check trial inclusion has been run already *****
  if ~isfield(Info,'Tlist')
     disp('Failed to identify trial inclusion criteria for use in Compute_DSI');
     disp('Make sure Trial_Inclusion is run to update Info before calling');
     return;
  end
  
  %****** compute the DSI ******
  sp = Exp.sp{Unit};
  %*******
  Info.StimInt = StimInt;
  Info.StimWin = StimWin;
  %****
  BigList = [Info.AttList ; Info.IgnAList ; Info.IgnBList];  %all valid trials
  OriInd = zeros(size(BigList));  % integer orientation values per trial
  OriDir = zeros(size(BigList));  % integer orientation values per trial
  OriSpk = zeros(size(BigList));  % spike counts per trial
  %***** for reference, what are the set of orientations
  OriVals = unique(Info.StimOri);  
  NVals = length(OriVals);
  %******** resort trials by order of orientation
  StimRast = [];  % keep raster of stim response
  %****** build a column of orienation values (int) and spike rates
  disp('Collecting spike counts per orientation for stimulus tuning ...');
  
  spike_times = sp.st;
  event_times = [];
  
  for k = 1:length(BigList)
      tr = BigList(k);
      trial = Info.Tlist(tr);
      %*************
      OriDir(k) = Info.StimOri(tr);
      kori = find(OriDir(k) == OriVals);
      OriInd(k) = kori;
      %****** get stim spike count on trial
      onstim = Info.StimOnList(tr);
      onstaspk = Exp.D{trial}.START_EPHYS + onstim + StimInt(1); % in secs
      onfinspk = Exp.D{trial}.START_EPHYS + onstim + StimInt(2);
      z = find( (sp.st >= onstaspk) & (sp.st < onfinspk) );
      if ~isempty(z)
           rate = (length(z)/(onfinspk-onstaspk));
      else
           rate = 0;
      end
      OriSpk(k) = rate;
      %****** store raster *****
      event_times = [event_times ; Exp.D{trial}.START_EPHYS + onstim];
      %***********
      staspk = Exp.D{trial}.START_EPHYS + onstim + StimWin(1); % in secs
      finspk = Exp.D{trial}.START_EPHYS + onstim + StimWin(2);
      z = find( (sp.st >= staspk) & (sp.st < finspk) );
      if ~isempty(z)
        difftime = (sp.st(z) - staspk );
        StimRast = [StimRast ; [(difftime+StimWin(1)) (ones(size(difftime))*k)]];
      end
      %*************************
  end 
  %**** store local info
 % data.spike_times = spike_times;
 % data.event_times = event_times;
 % save('psthdata','data');
  
  %****** compute the DSI with a sliding window over firing rate
  CountWin = 0.02;  % small counting bins for this, 20 ms width
  [tt,uu,su] = spikestats.CompCountSlider(StimRast,[],[],OriVals,OriInd,StimWin,CountWin,3);
  Info.DSI_tt = tt;
  Info.DSI_uu = uu;
  Info.DSI_su = su;
  %******** other computations then
  Info.StimRast = StimRast;
  Info.OriInd = OriDir;  % trial list of orientation per trial 
  VM = length(OriDir);   % must give PSTH number of trials to fill rast
  [tt,uu,su] = spikestats.CompPSTH(Info.StimRast,Info.StimWin,VM,Smooth);     
  Info.StimTT = tt;
  Info.StimUU = uu;   % figure rate before stimulus onset
  Info.StimSU = su;
  %******* compute the base firing rate in pre-stimulus period
   zz = find(tt < 0);  % get baseline firing as time before 0 ms
   if ~isempty(zz)
      Info.BaseMu = mean(uu(zz));
   else
      Info.BaseMu = NaN;
   end
   %***** determine if it meets inclusion test for excitatory vis response
   Info.VisCrit = 0;
   Info.SigTT = [];
   Info.ISigTT = [];
   if ~isnan(Info.BaseMu)
       vzz = find( (tt >= Info.StimInt(1)) & (tt < Info.StimInt(2)) );
       if ~isempty(vzz)
         testpts = length(vzz);
         sigo = norminv( 0.05 / testpts);  %Bonferonni correction
         sigt = [];
         isigt = [];
         for tk = 1:length(vzz)
            vmu = uu(vzz(tk));
            smu = su(vzz(tk));
            if ((vmu+(sigo*smu)) > Info.BaseMu)
               sigt = [sigt tt(vzz(tk))];
            end
            if ((vmu-(sigo*smu)) < Info.BaseMu)
               isigt = [isigt tt(vzz(tk))];
            end         
         end
         Info.SigTT = sigt;  % times vis response is significant
         Info.ISigTT = isigt;
         %**** test if more excitatory points than inhibitory
         % if (length(Info.SigTT) > length(Info.ISigTT))
         %    Info.VisCrit = 1;
         % end
         %***** next test, compute via spike counts
         vcounts = [];
         for zk = 1:length(Info.OriInd)
             if ~isempty( Info.StimRast )
               z = find( Info.StimRast(:,2) == zk );
               sptim = Info.StimRast(z,1); % times
               count = length( find( (sptim >= Info.StimInt(1)) & ...
                                   (sptim <  Info.StimInt(2)) ));
               vcounts = [vcounts ; (count / (Info.StimInt(2) - Info.StimInt(1)))];
             end
         end
         if isempty(vcounts)
             Info.VisMu = NaN;
             Info.VisSu = NaN;
             Info.VisCrit = 0;
         else
            Info.VisMu = mean(vcounts);
            Info.VisSu = std(vcounts) / sqrt(length(vcounts));
            %**** test if more excitatory points than inhibitory
            if ((Info.VisMu-(2*Info.VisSu)) > Info.BaseMu)
              % length(Info.SigTT) > length(Info.ISigTT))
              Info.VisCrit = 1;
            end
         end
         %********
       end     
   end
   %*************
  
  %****** compute the orientation tuning and DSI value
  %*****     and will implement as a Jacknife soon
  disp('Computing orientation tuning and DSI');
  Info.OriTune = spikestats.compute_orientation_dsi(OriVals,OriInd,OriSpk);
  BigList = [Info.AttList ; Info.IgnAList ; Info.IgnBList];  %all valid trials
  alist = (1:length(Info.AttList))';
  ilist = ((length(alist)+1):length(OriInd))';
  Info.AOriTune = spikestats.compute_orientation_dsi(OriVals,OriInd(alist),OriSpk(alist));
  Info.IOriTune = spikestats.compute_orientation_dsi(OriVals,OriInd(ilist),OriSpk(ilist));
  
  %****** take mean of att and ignored for pref (unbiased bet two)
  Wvec2 = (Info.AOriTune.Wvec/abs(Info.AOriTune.Wvec)) + ...
          (Info.IOriTune.Wvec/abs(Info.IOriTune.Wvec));
  Info.OriTune.Pref2 = spikestats.ResultantPref(OriVals,Wvec2);
  %*****
  
  %*** then run Jacknife on the DSI
  dsivals = [];
  angdev_vals = [];
  angstd_vals = [];
  NA = length(OriInd);
  for i = 1:JackN
      %****
      ia = 1 + (i-1)*floor(NA/JackN);
      ib = i * floor(NA/JackN);
      subset = [1:ia, ib:NA];
      %***
      tmp = spikestats.compute_orientation_dsi(OriVals,OriInd(subset),OriSpk(subset));
      dsivals = [dsivals tmp.DSI];
      angdev_vals = [angdev_vals tmp.AngDev];
      angstd_vals = [angstd_vals tmp.AngStd];
  end
  Info.OriTune.DSI_SEM = std(dsivals) * sqrt(JackN-1);
  Info.OriTune.ANGDEV_SEM = std(angdev_vals) * sqrt(JackN-1);
  Info.OriTune.ANGSTD_SEM = std(angstd_vals) * sqrt(JackN-1);
  %***** categorize neuron inclusion
  %Info.OriTune.LowDSI = Info.OriTune.DSI - (2 * Info.OriTune.DSI_SEM);
  %if (Info.OriTune.LowDSI >= DSI_threshold)
  if (Info.OriTune.DSI > DSI_threshold)   % allow more units, filter on fits
      Info.TunedUnit = 1;
  else
      Info.TunedUnit = 0;
  end
  %*********
  Info.VonMisesFit = 0;
  if ~isempty(Info.StimRast)
     disp('Fitting VonMises curve to ori tuning');    
     [Afit,Atune,~,~,~,~,~] = VonMises.FitWithVonMises(OriVals,JackN,OriDir,OriSpk,[]);
     if ~isempty(Afit)
       Info.OriTune.Afit = Afit;
       Info.OriTune.Atune = Atune;
       if (Afit.sem(4) < (THETAERROR*(pi/180)) )
          Info.VonMisesFit = 1;
       end
       %***** Determine VonMises PrefOri, call it VonPref (note should be
       %***** less noisy for neurons with low firing rates or low trials
       z = find( max(Atune.mu) == Atune.mu);
       Info.OriTune.VonPref = z(1);
       %*****************
     else
       Info.OriTune.Afit.mu = zeros(1,4);  % no values
       Info.OriTune.Afit.sem = zeros(1,4);  % no values
       Info.OriTune.Atune.mu = zeros(size(OriVals)); % zero tuning
       Info.OriTune.Atune.sem = ones(size(OriVals)); % zero tuning
     end
  else
     Info.OriTune.Afit.mu = zeros(1,4);  % no values
     Info.OriTune.Afit.sem = zeros(1,4);  % no values
     Info.OriTune.Atune.mu = zeros(size(OriVals)); % zero tuning
     Info.OriTune.Atune.sem = ones(size(OriVals)); % zero tuning
  end
  %****** Also consider trial numbers for inclusion in attention cases
  %****** and require some minimum number (48) per condition
  MIN_TRIALS = 48; %24;
  ATrials = length(Info.AttList);
  ITrials = length(Info.IgnAList)+length(Info.IgnBList);
  not_enough = 0;
  if (ATrials < MIN_TRIALS)
       not_enough = 1;
  end
  if (ITrials < MIN_TRIALS)
       not_enough = 1;
  end
  %******** Final decision to include neuron or not
  Info.Included = 0;
  if Info.VisCrit  % include even no sig DSI, throw out later from
                   % analyses based on Info.TunedUnit 
  %if  (Info.TunedUnit && Info.VisCrit) % && Info.VonMisesFit)
     if ~not_enough
         Info.Included = 1;
     end
  end
  %***********
  
return;


%******* below is the plot function in the same file
function plot_compute_dsi(Info,H)
   %******* This plots 1) firing rate amplitude as a function of motion direction
   %***                2) confidence intervals on that curve
   %***                3) computes DSI value with 95% confidence intervals
   
   subplot(H);
   OriTune = Info.OriTune;
   %*********
   px = []; py = [];
   px1 = []; py1 = [];
   px2 = []; py2 = [];
   rx = []; ry = [];
   for k = 1:length(OriTune.OriVals)
     ori = OriTune.OriVals(k)*(pi/180);
     vc = OriTune.Atune.mu(k) * complex(cos(ori),sin(ori));
     px = [px ; real(vc)];
     py = [py ; imag(vc)];
     %******
     vc = (OriTune.Atune.mu(k) + (2*OriTune.Atune.sem(k)) ) * complex(cos(ori),sin(ori));
     px1 = [px1 ; real(vc)];
     py1 = [py1 ; imag(vc)];
     %******
     vc = (OriTune.Atune.mu(k) - (2*OriTune.Atune.sem(k)) ) * complex(cos(ori),sin(ori));
     px2 = [px2 ; real(vc)];
     py2 = [py2 ; imag(vc)];
     %******
     vc = OriTune.Mu(k) * complex(cos(ori),sin(ori));
     rx = [rx ; real(vc)];
     ry = [ry ; imag(vc)];
  end
  px = [px ; px(1)];
  py = [py ; py(1)];
  px1 = [px1 ; px1(1)];
  py1 = [py1 ; py1(1)];
  px2 = [px2 ; px2(1)];
  py2 = [py2 ; py2(1)];
  rx = [rx ; rx(1)];
  ry = [ry ; ry(1)];  
  %**** draw the circle
  maxo = max(OriTune.Mu);
  %rmaxo = max(max(abs(rx)),max(abs(ry))) * 1.1;
  h = plot(rx,ry,'k.'); hold on;
  set(h,'Markersize',10);
  h = plot(px,py,'b-'); hold on;
  set(h,'Linewidth',2);
  plot(px1,py1,'b:'); hold on;
  plot(px2,py2,'b:'); hold on;
  h = plot([0,real(OriTune.Wvec * maxo)],[0,imag(OriTune.Wvec * maxo)],'r-');
  set(h,'Linewidth',2);
  h = plot(0,0,'ro');
  set(h,'Markersize',8);
  set(h,'Linewidth',2);
  axis tight;
  V = axis;
  rmaxo = max(abs(V)) * 1.1;
  axis([-rmaxo rmaxo -rmaxo rmaxo]);
  title(sprintf('DSI = %4.2f  conf( %4.2f to %4.2f )',OriTune.DSI,...
        (OriTune.DSI-(2*OriTune.DSI_SEM)),(OriTune.DSI+(2*OriTune.DSI_SEM))));    
  if (Info.Included)
     h = text((-rmaxo*0.15),(rmaxo*0.85),sprintf('INCLUDED ISO(%d)',Info.isolation(1)));
     set(h,'Fontsize',10);
     set(h,'Color',[1,0,0]);
  else
      if ~Info.VisCrit
          h = text((-rmaxo*0.15),(rmaxo*0.85),sprintf(' non-visual: ISO(%d)',Info.isolation(1)));
          set(h,'Fontsize',10);
          set(h,'Color',[0,0,1]);
      else
          if ~Info.VonMisesFit
             h = text((-rmaxo*0.15),(rmaxo*0.85),sprintf(' non-vonmise: ISO(%d)',Info.isolation(1)));
             set(h,'Fontsize',10);
             set(h,'Color',[0,1,0.5]);       
          else
             h = text((-rmaxo*0.15),(rmaxo*0.85),sprintf(' non-tuned: ISO(%d)',Info.isolation(1)));
             set(h,'Fontsize',10);
             set(h,'Color',[0,0.5,1]);
          end
      end
  end
  
return;

%******* below is the plot of orientation tuning
function plot_tuning(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   OriTune = Info.OriTune;
   spikeplot.TuningCurvePlot(OriTune,[0,0,0.6],'.',0,0);
   %****
   axis tight;
   V = axis;
   xlabel('Direction');
   ylabel('Rate');
   %*****
   if (V(4) > V(3))
     h3 = text((V(1)+0.2*(V(2)-V(1))),(V(3)+0.05*(V(4)-V(3))),...
             sprintf('Pref(%5.1f,sm:%5.1f)',(180/pi)*OriTune.Afit.mu(4),...
                                           (180/pi)*OriTune.Afit.sem(4)));
     set(h3,'Fontsize',10);
     set(h3,'Color',[0,0,1]);
     set(h3,'Fontweight','bold');
   end
   %*************
   title(sprintf('Name %s',Info.tagname));
   
return;


%******* below is the plot of stimulus locked PSTH
function plot_stimpsth(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   
   %****************
   VM = length(Info.OriInd);
   if isempty(Info.StimRast)   % no spikes from unit in the exp
       return;
   end
   %******************
   h = spikeplot.PlotTickRaster(Info.StimRast,Info.OriInd,0); hold on;
   set(h,'Linewidth',2);
   %*** superimpose a mean PSTH with error bars into same plot
   tt = Info.StimTT;
   uu = Info.StimUU;
   su = Info.StimSU;
   basemu = Info.BaseMu;
   viscrit = Info.VisCrit;
   %********* now plot results
   maxo = max( uu + (2*su));
   scal = (VM/maxo);
   h2 = plot(tt,uu*scal,'b-',tt,(uu+(2*su))*scal,'b-',tt,(uu-(2*su))*scal,'b-');
   %*** plot baseline rate line
   zz = find(tt < 0);  % get baseline firing as time before 0 ms
   if ~isnan(basemu)
      vmu = basemu*scal;
      plot([Info.StimWin(1),Info.StimWin(2)],[vmu,vmu],'b--');
      h3 = text((Info.StimWin(1)+0.10*(Info.StimWin(2)-Info.StimWin(1))),...
                 (vmu*1.0),sprintf('%4.1f hz',basemu));
      set(h3,'Fontsize',12);
      set(h3,'Color',[0,0,1]);
      set(h3,'Fontweight','bold');
   end
   %**************** super-impose DSI over time too *********
   dmax = max(Info.DSI_uu + (2*Info.DSI_su));
   if (dmax > 0)
       scal = VM;
       h3 = plot(Info.DSI_tt,scal*Info.DSI_uu,'g-');
       set(h3,'Linewidth',2);
       plot(Info.DSI_tt,scal*(Info.DSI_uu+(2*Info.DSI_su)),'g-');
       plot(Info.DSI_tt,scal*(Info.DSI_uu-(2*Info.DSI_su)),'g-');
   end
   %***************
   axis([Info.StimWin(1) Info.StimWin(2) 0 VM]);
   plot([0,0],[0,VM],'k-');
   plot([Info.StimInt(1),Info.StimInt(1)],[0,VM],'k--');
   plot([Info.StimInt(2),Info.StimInt(2)],[0,VM],'k--');
   vmm = 0.95 * VM;
   if ~isempty(Info.SigTT)
     plot(Info.SigTT,(vmm*ones(size(Info.SigTT))),'b.');
   end
   if ~isempty(Info.ISigTT)
     plot(Info.ISigTT,(vmm*ones(size(Info.ISigTT))),'c.');
   end
   xlabel('Time (secs)');
   ylabel('Trials');
   title(sprintf('Vis(%4.1f,%4.1f) z(%4.1f)',Info.VisMu,Info.VisSu,((Info.VisMu-Info.BaseMu)/Info.VisSu)));
   %******************
   
return;
   
function SInfo = Compute_Saccade_Lock_Forage_Motor(Info,Exp,Unit,H,TrialType,Interval)
%******* 
%******  function Info = Compute_DSI(Info,Exp,Unit,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          Unit - spike cluster to run computation on
%***          H - handle to plot window if provided, or if 1 create, else []
%***          TrialType - 1 use BackImage, 2 use Forage (noisetype = 1),
%***                        3 use Forage (noisetype = 3), 
%***                        4 use Forage (noisetype = 6)
%**           Interval - use this time after saccade onset as a window
%**                       to compute the motor tuning, and if multiple
%**                       columns then do each interval
%***
%*** Outputs: Info - it computes the DSI from all trial types,
%***                   and puts confidence intervals on it

  SInfo = [];  % return saccade Info object
  
  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh = figure('position',[100 100 800 800]); 
       HA = subplot('position',[0.1 0.35 0.5 0.6]);
       HB = subplot('position',[0.1 0.10 0.5 0.2]);
       HC = subplot('position',[0.65 0.70 0.25 0.25]);
       HD = subplot('position',[0.65 0.40 0.25 0.25]);
       HE = subplot('position',[0.65 0.10 0.25 0.25]);
    else
        if isempty(H)
          return;
        else
          if (length(H) ~= 5)
              disp('Error with subplot pass to Compute_Saccade_Lock, requires two subplots');
              return;
          else
              HA = H{1};
              HB = H{2};
              HC = H{3};
              HD = H{4};
              HE = H{5};
          end
        end
    end
    %*******
    plot_sacrast(Info,HA);
    plot_sacpsth(Info,HB);
    plot_sacdist(Info,HC);
    plot_sactune(Info,HD,1);  % tuning first interval 
    plot_sactune(Info,HE,2);  % tuning second interval
    
    %********
    if isfloat(H)  && (H == 2)
       z = getframe(hh);  % current figure
       uname = [Info.pathplot,filesep,sprintf('SAC%d_',TrialType),Info.tagname,'.png'];
       disp(sprintf('Storing image of SAC graph at %s',uname)); 
       imwrite(z.cdata,uname); % store image for later review
       close(hh);
    end   
    return;
  end
  
  %****** parameters for computing DSI and neuron inclusion
  JackN = 10;             % Jacknife to get confidence intervals
  StimWin = [-0.1,0.2];   % plot stim locked firing rate PSTH 
  Smooth = 5;             % Gaussian sig for smooth PSTH plots
  
  %****** use trial inclusion to build list of trials  *****
  TList = [];
  for fr=1:size(Exp.D,1)
    if (TrialType == 2)  
      if strcmp(Exp.D{fr}.PR.name,'ForageProceduralNoise') 
        if (Exp.D{fr}.PR.noisetype == 1)
           TList = [TList ; fr];             
        end
      end
    elseif (TrialType == 3)
        if strcmp(Exp.D{fr}.PR.name,'Forage') 
          if (Exp.D{fr}.PR.noisetype == 3)
           TList = [TList ; fr];             
          end  
        end
    elseif (TrialType == 4)
        if strcmp(Exp.D{fr}.PR.name,'Forage') 
          if (Exp.D{fr}.PR.noisetype == 6)
           TList = [TList ; fr];             
          end  
        end
    else
        if strcmp(Exp.D{fr}.PR.name,'BackImage') 
           TList = [TList ; fr];             
        end        
        
    end
  end
  %******* then go to collect all the saccade events
  EventList = [];  % list of times in ephys time for time-locking, and 
                   % second column is RT to next saccade
  %*** store saccade vector per event, and let's bin them too
  SacVecList = [];  % [x,y] vector (column 1+2), bin - integer last column
  N_SacDir = 12;  % bins to code direction
  N_SacEcc = 3;  % bins to code eccentricity, power of 2, thus <2, < 4, <8, <16
  %*******
  for tr = 1:length(TList)
      fr = TList(tr);  % trial number
      estart = Exp.D{fr}.START_EPHYS;
      % MatStart = Exp.D{fr}.eyeData(6,1);  
      for kt = 1:size(Exp.D{fr}.slist,1)  % list of flagged saccades       
         if (kt < size(Exp.D{fr}.slist,1) )
            % RT = (Exp.D{fr}.slist(kt+1,1) - Exp.D{fr}.slist(kt,1));
            RT = (Exp.D{fr}.slist(kt+1,2) - Exp.D{fr}.slist(kt,2));
         else
            RT = 1;  % else set as long to next saccade, 1 second
         end
         etime = estart + Exp.D{fr}.slist(kt,1);  % saccade onset
          % etime = estart + Exp.D{fr}.slist(kt,2);
         EventList = [EventList ; [etime,RT]];
         
         %******** saccade information
         sta = Exp.D{fr}.slist(kt,4); % integer, sac start
         eta = Exp.D{fr}.slist(kt,5); % integer, sac end
         dx = Exp.D{fr}.eyeSmo(eta,2) - Exp.D{fr}.eyeSmo(sta,2);
         dy = Exp.D{fr}.eyeSmo(eta,3) - Exp.D{fr}.eyeSmo(sta,3);
         mag = norm([dx,dy]);
         ang = (angle(complex(dx,dy))+pi)*(180/pi);
         bin = ( (floor(log2(max(mag,1)))) * N_SacDir);
         bin = bin + 1 + floor(ang/(360/N_SacDir));
         SacVecList = [SacVecList ; [dx,dy,mag,ang,bin]];
         %***************
      end
  end
  
  %******** no events found, then return
  if isempty(EventList) || isempty(SacVecList)
      return;
  end
  
  %****** compute the raster and PSTH ******
  sp = Exp.sp{Unit};
  %*******
  SInfo.EventList = EventList;
  SInfo.StimWin = StimWin;
  SInfo.TrialType = TrialType;
  SInfo.SacVecList = SacVecList;
  SInfo.N_SacDir = N_SacDir;
  SInfo.N_SacEcc = N_SacEcc;
  %****
  %******** resort trials by order of orientation
  StimRast = [];  % keep raster of stim response
  %****** build a column of orienation values (int) and spike rates
  disp('Collecting spike counts per saccade  ...');
  
  spike_times = sp.st;
  event_times = [];
  for k = 1:length(EventList)
      OriInd(k) = EventList(k,2);  % use RT time to sort tirals
      %****** get stim spike count on trial
      onstim = EventList(k,1); 
      staspk = onstim + StimWin(1); % in secs
      finspk = onstim + StimWin(2);
      z = find( (sp.st >= staspk) & (sp.st < finspk) );
      if ~isempty(z)
        difftime = (sp.st(z) - staspk );
        StimRast = [StimRast ; [(difftime+StimWin(1)) (ones(size(difftime))*k)]];
      end
      %*************************
  end 
  if isempty(StimRast)
      SInfo = [];
      return;
  end
  
  %******** other computations then
  SInfo.StimRast = StimRast;
  SInfo.OriInd = EventList(:,2);  % trial list of RTs per event 
  VM = length(SInfo.OriInd);   % must give PSTH number of trials to fill rast
  [tt,uu,su] = spikestats.CompPSTH(SInfo.StimRast,SInfo.StimWin,VM,Smooth);     
  SInfo.StimTT = tt;
  SInfo.StimUU = uu;   % figure rate before stimulus onset
  SInfo.StimSU = su;
  %******* compute the base firing rate in pre-stimulus period
  zz = find(tt < 0);  % get baseline firing as time before 0 ms
  if ~isempty(zz)
      SInfo.BaseMu = mean(uu(zz));
  else
      SInfo.BaseMu = NaN;
  end
  %*************
  
  %******** now run Sac Events, get spike counts
  NI = size(Interval,1);
  if NI > 0
    SacTune = cell(1,NI);
    SacUU = cell(1,NI);
    SacSU = cell(1,NI);
    SacNN = cell(1,NI);
    for nk = 1:NI
      SacTune{nk} = cell(N_SacEcc,N_SacDir);
      SacUU{nk} = NaN(N_SacEcc,N_SacDir);
      SacSU{nk} = NaN(N_SacEcc,N_SacDir);
      SacNN{nk} = NaN(N_SacEcc,N_SacDir);
    end
    for k = 1:length(EventList)
      sbin = SacVecList(k,5);  % bin number
      secc = 1+floor((sbin-1)/N_SacDir);
      if (secc > N_SacEcc)
          continue;  % outside bounds
      else
          sdir = 1+mod((sbin-1),N_SacDir);
          %****** get spike count, 1st time bin
          for nk = 1:NI 
            onstim = EventList(k,1); 
            staspk = onstim + Interval(nk,1); % in secs
            finspk = onstim + Interval(nk,2);
            z = find( (sp.st >= staspk) & (sp.st < finspk) );
            pcount = (length(z)/(Interval(nk,2)-Interval(nk,1))); 
            SacTune{nk}{secc,sdir} = [SacTune{nk}{secc,sdir} ; pcount];
          end
          %**************
      end
    end
  end
  %********* compute mean and sem on the saccade tuning curves
  for nk = 1:NI
    for i = 1:N_SacEcc
      for j = 1:N_SacDir
          %******* summarize for peak values
          SacNN{nk}(i,j) = length(SacTune{nk}{i,j});
          if ~isempty(SacTune{nk}{i,j})
             SacUU{nk}(i,j) = nanmean(SacTune{nk}{i,j});
             SacSU{nk}(i,j) = nanstd(SacTune{nk}{i,j})/sqrt(SacNN{nk}(i,j));
          end
      end
    end
  end
  %***
  SInfo.SacTune = SacTune;
  SInfo.SacUU = SacUU;
  SInfo.SacSU = SacSU;
  SInfo.SacNN = SacNN;
  SInfo.NI = NI;
  SInfo.Interval = Interval;
  %*********** Now, we could condense into circular tuning curve too here? 
  SACUR = cell(NI,2);
  SACDSI = NaN(NI,1);
  SACDSISTD = NaN(NI,1);
  SACDSISHUF = NaN(NI,1);
  SACDSIPVAL = NaN(NI,1);
  %**********
  for nk = 1:NI
          %******* get the raw DSI metric with boot-strap on it
          zOriVals = [];
          zOriInd = [];
          zOriSpk = [];
          for ii = 1:SInfo.N_SacDir
              zOriVals(ii) = (ii-1)*(360/SInfo.N_SacDir);
              for jj = 1:SInfo.N_SacEcc
                  vals = SInfo.SacTune{nk}{jj,ii};
                  if ~isempty(vals)
                    zOriSpk = [zOriSpk ; vals];
                    zOriInd = [zOriInd ; (zOriVals(ii)*ones(size(vals)))];
                  end
              end
          end
          rr = randperm(length(zOriInd));
          OriInd = zOriInd(rr);
          OriSpk = zOriSpk(rr);
          OriVals = zOriVals; % unique(OriInd)';
          %********
          BootSamp = 1000;
          dsivals = [];
          dsishuf = [];
          NA = length(OriInd);
          disp('Running bootstrap for saccade tuning ...');
          for i = 1:BootSamp
              %****
              rp = randi(length(OriInd),length(OriInd),1);
              rr = randperm(length(OriInd)); % a reshuffle of orientation indices
              %*********
              tmp = spikestats.compute_orientation_dsi2(OriVals,OriInd(rp),OriSpk(rp));
              dsivals = [dsivals tmp.DSI];
              %****** now get the shuffle corrector
              tmp = spikestats.compute_orientation_dsi2(OriVals,OriInd(rr),OriSpk(rp));
              dsishuf = [dsishuf tmp.DSI];    
          end
          disp('... bootstrap finished');
          SACDSI(nk,1) = mean(dsivals);
          SACDSISTD(nk,1) = std(dsivals);
          SACDSISHUF(nk,1) = mean(dsishuf);
          zz = find( dsishuf > SACDSI(nk,1)); 
          SACDSIPVAL(nk,1) = length(zz)/BootSamp;  % do this way, not normal distribution
          %*** compute column averages and SEM for plots
          sacu = [];
          sacsem = [];
          if (1)
            for sk = 1:SInfo.N_SacDir
                zz = find( OriInd == OriVals(sk));
                uu = nanmean(OriSpk(zz));
                su = nanstd(OriSpk(zz))/sqrt(length(OriSpk(zz)));
                sacu = [sacu ; uu];
                sacsem = [sacsem ; su];
            end
          end
          SACUR{nk,1} = [SACUR{nk,1} ; sacu];
          SACUR{nk,2} = [SACUR{nk,2} ; sacsem];
          %******************
          maxo = max(sacu);  % avg over SF, take peak
          mino = min(sacu);
          AImot = (maxo-mino)/(maxo+mino);  % range 0 to 1
          %******** look at example units
          if (0)  % for debugging
             H = figure(11);
             set(H,'Position',[100 550 350 350]);
             xx = (1:SInfo.N_SacDir)*(360/SInfo.N_SacDir);
             plot(xx,sacu,'k-','Linewidth',2); hold on;
             plot(xx,(sacu+(2*sacsem)),'k-');
             plot(xx,(sacu-(2*sacsem)),'k-');
             plot(xx,(mean(sacu)*ones(size(sacu))),'k:');
             xlabel('Saccade Direction');
             ylabel('Rate (sp/s)');
             title(sprintf('AI %5.3f DSI %5.3f(%6.4f)',AImot,SACDSI(nk,1),SACDSIPVAL(nk,1)));
             input('check unit');
             close all;
          end
  end % for nk = 1:NI
  %************
  SInfo.SACDSI = SACDSI;
  SInfo.SACDSISTD = SACDSISTD;
  SInfo.SACDSISHUF = SACDSISHUF;
  SInfo.SACDSIPVAL = SACDSIPVAL;
  SInfo.SACUR = SACUR;
  %************  
return;

%******* below is the plot of stimulus locked PSTH
function plot_sacrast(Info,H)
   %******* This plots 1) firing rate as function of saccade locking,
   %sorted by the RT (RT stored in OriInd)
   subplot(H);
   %****************
   if isempty(Info)
       return;
   end
   if isempty(Info.StimRast)   % no spikes from unit in the exp
       return;
   end
   VM = length(Info.OriInd);
   %******************
   h = spikeplot.PlotTickRaster(Info.StimRast,Info.OriInd,0); hold on;
   set(h,'Linewidth',2);
   %*** superimpose a mean PSTH with error bars into same plot
   tt = Info.StimTT;
   uu = Info.StimUU;
   su = Info.StimSU;
   basemu = Info.BaseMu;
   %***************
   axis([Info.StimWin(1) Info.StimWin(2) 0 VM]);
   plot([0,0],[0,VM],'k-');
   xlabel('Time (secs)');
   ylabel('Trials');
   if (Info.TrialType == 1)
       title('Background Images');
   else
       if (Info.TrialType == 2)
           title('Grating Stimuli');
       else
          if (Info.TrialType == 3)
             title('Blank Screen');
          else
             title('Full-field Dots')       
          end
       end
   end
   %******************
   
return;

function plot_sacpsth(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   %*******
   if isempty(Info)
       return;
   end
   if isempty(Info.StimRast)   % no spikes from unit in the exp
       return;
   end
   %*** superimpose a mean PSTH with error bars into same plot
   tt = Info.StimTT;
   uu = Info.StimUU;
   su = Info.StimSU;
   basemu = Info.BaseMu;
   %********* now plot results
   VM = max( uu + (2*su))*1.1;
   h2 = plot(tt,uu,'b-',tt,(uu+(2*su)),'b-',tt,(uu-(2*su)),'b-'); hold on;
   %*** plot baseline rate line
   zz = find(tt < 0);  % get baseline firing as time before 0 ms
   if ~isnan(basemu)
      vmu = basemu;
      plot([Info.StimWin(1),Info.StimWin(2)],[vmu,vmu],'b--');
      h3 = text((Info.StimWin(1)+0.10*(Info.StimWin(2)-Info.StimWin(1))),...
                 (vmu*1.0),sprintf('%4.1f hz',basemu));
      set(h3,'Fontsize',12);
      set(h3,'Color',[0,0,1]);
      set(h3,'Fontweight','bold');
   end
   %***************
   axis([Info.StimWin(1) Info.StimWin(2) 0 VM]);
   plot([0,0],[0,VM],'k-');
   xlabel('Time (secs)');
   ylabel('Rate (sp/s)');
   %**********
   
return;

%********** Plot motor tuning at peak and dip, z-scored rates
function plot_sacdist(Info,H)
   subplot(H);
   plot(Info.SacVecList(:,1),Info.SacVecList(:,2),'k.','Markersize',2); hold on;
   axis([-16 16 -16 16]);
return;

function plot_sactune(Info,H,Inter)  
   %**********
   subplot(H);
   upp = nanmean(nanmean(Info.SacUU{Inter}));
   zpeak = (Info.SacUU{Inter} - upp) ./ Info.SacSU{Inter};  % zscore effect
   imagesc(zpeak,[-8 8]); hold on;
   colorbar;
   title(sprintf('Interval %d',Inter));
return;

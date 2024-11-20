function SInfo = Compute_Saccade_Lock(Info,Exp,Unit,H,TrialType)
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
%***
%*** Outputs: Info - it computes the DSI from all trial types,
%***                   and puts confidence intervals on it

  SInfo = [];  % return saccade Info object
  
  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh = figure('position',[100 100 400 600]); 
       HA = subplot('position',[0.1 0.35 0.8 0.6]);
       HB = subplot('position',[0.1 0.10 0.8 0.2]);
    else
        if isempty(H)
          return;
        else
          if (length(H) ~= 2)
              disp('Error with subplot pass to Compute_Saccade_Lock, requires two subplots');
              return;
          else
              HA = H{1};
              HB = H{2};
          end
        end
    end
    %*******
    plot_sacrast(Info,HA);
    plot_sacpsth(Info,HB);
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
         etime = estart + Exp.D{fr}.slist(kt,1); %saccade onset
          % etime = estart + Exp.D{fr}.slist(kt,2);
         EventList = [EventList ; [etime,RT]];
      end
  end
  
  
  %****** compute the raster and PSTH ******
  sp = Exp.sp{Unit};
  %*******
  SInfo.EventList = EventList;
  SInfo.StimWin = StimWin;
  SInfo.TrialType = TrialType;
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
return;

%******* below is the plot of stimulus locked PSTH
function plot_sacrast(Info,H)
   %******* This plots 1) firing rate as function of saccade locking,
   %sorted by the RT (RT stored in OriInd)
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
   %***************
   axis([Info.StimWin(1) Info.StimWin(2) 0 VM]);
   plot([0,0],[0,VM],'k-');
   xlabel('Time (secs)');
   ylabel('Trials');
   if (Info.TrialType == 1)
       title('Background Images');
   else
       title('Grating Stimuli');
   end
   %******************
   
return;

function plot_sacpsth(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
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

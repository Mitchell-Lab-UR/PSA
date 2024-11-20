function Info = Compute_TimeLock_MatchRT(Info,Exp,Unit,H,FieldName,Event,WinLock,TimeLock,UseOri,MatchNum)
%******* 
%******  function Info = Compute_StimLock(Info,Exp,Unit,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          Unit - spike cluster to run computation on
%***          H - handle to plot window if provided, or if 1 create, else []
%***
%***          FieldName - create a new field with variable in it of name
%***          Event - if 1, use stim onset to time lock, 
%***                  if 2, use saccade onset
%***                  if 3, use saccade offset
%***          WinLock - window of analyses
%***          TimeLock - window of spike counting
%***          UseOri - use orientation for rast sorting
%***          MatchNum - resamples from distribution to averaging
%***
%*** Outputs: Info - it computes the stim locked PSTH and rasters
%***               - it also fits a VonMises to stim-locked attention conditions


  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh = figure('position',[200 70 800 800]); 
       HA = subplot('position',[0.075 0.60 0.375 0.35]);  % PSTH
       HB = subplot('position',[0.075 0.10 0.375 0.35]);  % Tuning Curve
       HC = subplot('position',[0.575 0.10 0.375 0.85]);  % Raster
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
    plot_timelock_psth(Info,FieldName,HA);
    plot_timelock_tuning(Info,FieldName,HB);
    plot_timelock_rast(Info,FieldName,HC);  % need to figure out how I want this
    %********
    if isfloat(H)  && (H == 2)
       z = getframe(hh);  % current figure
       uname = [Info.pathplot,filesep,FieldName,'_TimeLock_',Info.tagname,'.png'];
       disp(sprintf('Storing image of TimeLock graph at %s',uname)); 
       imwrite(z.cdata,uname); % store image for later review
       close(hh);
    end   
    return;
  end
  
  %****** parameters for computing stimlock response 
  tinfo = []; % new struct for this local 'FieldName'
  tinfo.JackN = 10;           % Jacknife to get confidence intervals
  tinfo.LockInt = TimeLock;   % Interval to compute motion tuning  
  tinfo.LockWin = WinLock;    % plot stim locked firing rate PSTH 
  tinfo.Smooth = 5;           % Gaussian sig for smooth PSTH plots
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
  OriVals = unique(Info.StimOri);  
  NVals = length(OriVals);
  
  % now passed in as an argument *****
  %MatchNum = 10;   % run the analysis 20 times, each time matching the
                   % attended and ignored RT distributions, and then
                   % pool all the results 
  minfo = cell(1,MatchNum);
  for mm = 1:(MatchNum+1)
      
    %****
    AttList = Info.AttList;
    IgnList = [ Info.IgnAList ; Info.IgnBList ];
    if (mm > 1)
        RTList = Info.SacOnList - Info.StimOnList;
        BigList = spikestats.Match_Distribution_Lists(AttList,IgnList,RTList,0.01);
    else
        BigList = cell(1,2);
        BigList{1} = AttList;
        BigList{2} = IgnList;
    end
    %*****
    Rast = cell(1,2);
    binfo = [];   % struct to store cumulative results per match run
    %******
    for zk = 1:2 
      OriInd{zk} = zeros(size(BigList{zk}));  % integer orientation values per trial
      OriDir{zk} = zeros(size(BigList{zk}));  % integer orientation values per trial
      OriSpk{zk} = zeros(size(BigList{zk}));  % spike counts per trial
      OriRT{zk}  = zeros(size(BigList{zk}));  % list of reaction times
      %***** for reference, what are the set of orientations
      %******** resort trials by order of orientation
      %****** build a column of orienation values (int) and spike rates
      disp('Collecting spike counts per orientation for stimulus tuning ...');
      for k = 1:length(BigList{zk})
        tr = BigList{zk}(k);
        trial = Info.Tlist(tr);
        %*************
        OriDir{zk}(k) = Info.StimOri(tr);
        kori = find(OriDir{zk}(k) == OriVals);
        OriInd{zk}(k) = kori;
        %****** get stim spike count on trial
        if (Event == 1)
         tlock = Info.StimOnList(tr);
        end
        if (Event == 2);
         tlock = Info.SacOnList(tr);
        end
        if (Event == 3)
         tlock = Info.SacOffList(tr);
        end
        %***********
        OriRT{zk}(k) = (Info.SacOnList(tr) - Info.StimOnList(tr));
        onstaspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockInt(1); % in secs
        onfinspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockInt(2);
        z = find( (sp.st >= onstaspk) & (sp.st < onfinspk) );
        if ~isempty(z)
           rate = (length(z)/(onfinspk-onstaspk));
        else
           rate = 0;
        end
        OriSpk{zk}(k) = rate;   
        %****** store raster *****
        staspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockWin(1); % in secs
        finspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockWin(2);
        z = find( (sp.st >= staspk) & (sp.st < finspk) );
        if ~isempty(z)
         difftime = (sp.st(z) - staspk );
         Rast{zk} = [Rast{zk} ; [(difftime+tinfo.LockWin(1)) (ones(size(difftime))*k)]];
        end
      end
      %*************************
      if (mm == 1)
         %***** save full descripts for plotting rasters etc... (makes no sense to average rasters) 
         tinfo.Event = Event;            % record what you time lock on
         tinfo.UseOri = UseOri;          % use orientation to raster index
         tinfo.Rast{zk} = Rast{zk};
         tinfo.OriInd{zk} = OriDir{zk};  % trial list of orientation per trial (could do RT instead)
         tinfo.OriIndRT{zk} = OriRT{zk}; % need to figure this out
      end
      VM = length(OriDir{zk});   % must give PSTH number of trials to fill rast
      [tt,uu,su] = spikestats.CompPSTH(Rast{zk},tinfo.LockWin,VM,tinfo.Smooth);     
      binfo.LockTT{zk} = tt;
      binfo.LockUU{zk} = uu;
      binfo.LockSU{zk} = su;
      %****** compute the orientation tuning and DSI value
      %*****     and will implement as a Jacknife soon
      disp(sprintf('Computing orientation tuning and DSI: Attention Set(%d)',zk));
      otune = spikestats.compute_orientation_dsi(OriVals,OriInd{zk},OriSpk{zk});
      if (mm == 1)
          tinfo.OriTune{zk} = otune;
      else
          binfo.OriTune{zk}.Mu = otune.Mu;   % select some variables to average RT matched (not DSI)
          binfo.OriTune{zk}.Num = otune.Num;
          binfo.OriTune{zk}.Sem = otune.Sem;
      end
      if (mm == 1)
        %*** then run Jacknife on the DSI
        dsivals = [];
        angdev_vals = [];
        angstd_vals = [];
        NA = length(OriInd{zk});
        for i = 1:tinfo.JackN
          %****
          ia = 1 + (i-1)*floor(NA/tinfo.JackN);
          ib = i * floor(NA/tinfo.JackN);
          subset = [1:ia, ib:NA];
          %***
          tmp = spikestats.compute_orientation_dsi(OriVals,OriInd{zk}(subset),OriSpk{zk}(subset));
          dsivals = [dsivals tmp.DSI];
          angdev_vals = [angdev_vals tmp.AngDev];
          angstd_vals = [angstd_vals tmp.AngStd];
        end
        tinfo.OriTune{zk}.DSI_SEM = std(dsivals) * sqrt(tinfo.JackN-1);
        tinfo.OriTune{zk}.ANGDEV_SEM = std(angdev_vals) * sqrt(tinfo.JackN-1);
        tinfo.OriTune{zk}.ANGSTD_SEM = std(angstd_vals) * sqrt(tinfo.JackN-1);
      end
      %*************
      disp(sprintf('Fitting Unconstrained VonMises curve to ori tuning: Attention Set(%d)',zk));    
      % ofit = Info.OriTune.Afit.mu;
      % ofit = [ofit NaN];  % if an extra point with NaN, use start point, but do not constrain ori pref
      [AAfit,AAtune,~,~,~,~,~] = VonMises.FitWithVonMises(OriVals,tinfo.JackN,OriDir{zk},OriSpk{zk},[]);
      if ~isempty(AAfit)
        binfo.OriTune{zk}.AAfit = AAfit;
        binfo.OriTune{zk}.AAtune = AAtune;
      else
        binfo.OriTune{zk}.AAfit.mu = zeros(1,4);  % no values
        binfo.OriTune{zk}.AAfit.sem = zeros(1,4);  % no values
        binfo.OriTune{zk}.AAtune.mu = zeros(size(OriVals)); % zero tuning
        binfo.OriTune{zk}.AAtune.sem = ones(size(OriVals)); % zero tuning
      end
      %***** Then perform a fit with orientation constrained to pref ori    
      disp(sprintf('Fitting Constrained VonMises curve to ori tuning: Attention Set(%d)',zk));    
      Ofit = Info.OriTune.Afit.mu;  %constrain fit to start from DSI fit curve and fix pref ori, else do []
      [Afit,Atune,~,~,~,~,~] = VonMises.FitWithVonMises(OriVals,tinfo.JackN,OriDir{zk},OriSpk{zk},Ofit);
      if ~isempty(Afit)
        binfo.OriTune{zk}.Afit = Afit;
        binfo.OriTune{zk}.Atune = Atune;
      else
        binfo.OriTune{zk}.Afit.mu = zeros(1,4);  % no values
        binfo.OriTune{zk}.Afit.sem = zeros(1,4);  % no values
        binfo.OriTune{zk}.Atune.mu = zeros(size(OriVals)); % zero tuning
        binfo.OriTune{zk}.Atune.sem = ones(size(OriVals)); % zero tuning
      end
    end   % loop over zk conditions
    
    %****** store binfo to minfo, or otherwise set into tinfo as separate
    %****** struct containing the full fits with unmatched distribution
    %****** as a reference (how much does matching change results?)
    if (mm > 1)
          minfo{mm-1} = binfo;
    else
          tinfo.NonMatch = binfo;
    end
    %*************
  end
  %******** AVERAGE OVER RT matching
  tinfo = spikestats.Average_Matched_Fields_TimeLock(tinfo,minfo);  % average all subfields minfo, store results into tinfo
  %**********************************
  
  %*********** save it to a new field *******
  Info = setfield(Info,FieldName,tinfo);
  %*********
  
return;


%******* below is the plot of orientation tuning
function plot_timelock_tuning(Info,Field,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   tinfo = getfield(Info,Field);
   %****** need to show three plots potentially, but for now just two 
   for zk = 1:2    %could do 1:3 if want to see less sampled ignored loc
     if (zk == 1)
         colo = [1,0,0]; % att in red
         symb = '.';
         xoff = -2;
     else
         colo = [0,0,1]; % ign1 in blue
         symb = '.';
         xoff = 2;
     end
     OriTune = tinfo.OriTune{zk};  % code below is the same but for color and 
                                      % and input tuning, (maybe sep func)
     spikeplot.TuningCurvePlot(OriTune,colo,symb,xoff,1);
   end
   %*********
   axis tight;
   V = axis;
   xlabel('Direction (degs)','FontSize',15);
   ylabel('Rate (sp/s)','FontSize',15);
   title(sprintf('Unit Iso %d',Info.isolation(1)),'FontSize',18);
   % ax = gca;
   % ax.XAxis.FontSize = 15;
   % ax.YAxis.FontSize = 15;
   %**********
   
return;

%******* below is the plot of stimulus locked PSTH
function plot_timelock_psth(Info,Field,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   subplot(H);
   tinfo = getfield(Info,Field);
   spikeplot.plot_attention_traces(tinfo.LockTT,tinfo.LockUU,tinfo.LockSU,...
                                      tinfo.LockWin,tinfo.LockInt,tinfo.Event);
   ylabel('Rate (sp/s)','FontSize',15');
   title(sprintf('Name %s',Info.tagname),'FontSize',18);
   % ax = gca;
   % ax.XAxis.FontSize = 15;
   % ax.YAxis.FontSize = 15;
   %*************
return;

%******* below is the plot of stimulus locked PSTH
function plot_timelock_rast(Info,Field,H)
   %******* This plots raster plot across att conditions, stim locked
   subplot(H);
   tinfo = getfield(Info,Field);
   %******
   for zk = 1:2    %could do 1:3 if want to see less sampled ignored loc
     if (zk == 1)
         colo = [1,0,0]; % att in red
         voffset = 0;
     else
         colo = [0,0,1]; % ign1 in blue
         voffset = length(tinfo.OriInd{1});  % offset from first
     end
     %***
     if (tinfo.UseOri == 1)
       % arrange trials by orientation
       h2 = spikeplot.PlotTickRaster(tinfo.Rast{zk},tinfo.OriInd{zk},voffset); hold on;
     else
       % arrange trials by reaction time
       h2 = spikeplot.PlotTickRaster(tinfo.Rast{zk},tinfo.OriIndRT{zk},voffset); hold on;
     end
     set(h2,'Linewidth',2);
     set(h2,'Color',colo);
   end
   axis tight;
   V = axis;
   VM = V(4);
   axis([tinfo.LockWin(1) tinfo.LockWin(2) 0 VM]);
   plot([0,0],[1,VM],'k-');
   if (1) % use gray fill to mark time window
      aa = [tinfo.LockInt(1),tinfo.LockInt(1),tinfo.LockInt(2),tinfo.LockInt(2)];
      bb = [1,VM,VM,1];
      fill(aa,bb,[0.5,0.5,0.5],'FaceAlpha',0.3,'Linestyle','none');
   else
      plot([tinfo.LockInt(1),tinfo.LockInt(1)],[1,VM],'k--');
      plot([tinfo.LockInt(2),tinfo.LockInt(2)],[1,VM],'k--');
   end
   xlabel('Time (secs)','FontSize',15');
   if (tinfo.UseOri == 1)
      ylabel('Trials (sorted by Orientation)','FontSize',15');
   else
      ylabel('Trials (sorted by RT)','FontSize',15');       
   end
   
   % ax = gca;
   % ax.XAxis.FontSize = 15;
   % ax.YAxis.FontSize = 15;

   %******************
return;

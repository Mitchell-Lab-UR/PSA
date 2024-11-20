function Info = Compute_Discrim(Info,Exp,Unit,H,FieldName,Event,WinLock,TimeLocks)
%******* 
%******  function Info = Compute_Discrim(Info,Exp,Unit,H)
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
%***          TimeLocks - window of spike counting, cell array, each two
%***          RTpass - if not empty, sets min and max RT to include 
%***
%*** Outputs: Info - it computes the stim locked PSTH and rasters
%***               - it also fits a VonMises to stim-locked attention conditions


  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh = figure('position',[200 70 800 800]); 
       HA = subplot('position',[0.075 0.60 0.375 0.35]);  % PSTH
       HB = subplot('position',[0.075 0.10 0.375 0.35]);  % AUC
       HC = subplot('position',[0.575 0.60 0.375 0.35]);  % Fano
       HD = subplot('position',[0.575 0.10 0.375 0.35]);  % R Ori
    else
        if isempty(H)
          return;
        else
          if (length(H) ~= 2)
              disp('Error with subplot pass to Compute_Discrim, requires two subplots');
              return;
          else
              HA = H{1};
              HB = H{2};
              HC = H{3};
              HD = H{4};
          end
        end
    end
    %******* compute several stats over a sliding window
    plot_timelock_psth(Info,FieldName,HA);  % plot mean rate again
    plot_timelock_auc(Info,FieldName, HB);  % AUC over time
    % plot_timelock_fano(Info,FieldName, HC);  % plot FF over time
    plot_timelock_dsi(Info,FieldName, HD);  % plot ori std over time
    %********
    if isfloat(H)  && (H == 2)
       z = getframe(hh);  % current figure
       uname = [Info.pathplot,filesep,FieldName,'_Discrim_',Info.tagname,'.png'];
       disp(sprintf('Storing image of TimeLock graph at %s',uname)); 
       imwrite(z.cdata,uname); % store image for later review
       close(hh);
    end   
    return;
  end
  
  %****** parameters for computing stimlock response 
  tinfo = []; % new struct for this local 'FieldName'
  tinfo.JackN = 10;              % Jacknife to get confidence intervals
  tinfo.LockInt = TimeLocks;      % Interval to compute motion tuning  
  tinfo.LockWin = WinLock;       % plot stim locked firing rate PSTH 
  tinfo.Smooth = 5;              % Gaussian sig for smooth PSTH plots
  tinfo.CountWin = 0.05;         % Sliding window (secs) on AUC, Fano, OriStats
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
  
  %******* Determine the pref ori in computing AUC 
  % POri = Info.OriTune.Pref;
  % POri = Info.OriTune.VonPref;
  POri = Info.OriTune.Pref2;  % unbiased votes from att and ignored
  %*** need to derive what are three most and least pref directions
  PrefOri = [];
  NPrefOri = [];
  %*** for now just guess
  PrefOri = [POri-1,POri,POri+1];  % best three 
  z = find(PrefOri > length(OriVals));
  if ~isempty(z)
     PrefOri(z) = PrefOri(z) - length(OriVals);
  end
  z = find(PrefOri <= 0);
  if ~isempty(z)
     PrefOri(z) = PrefOri(z) + length(OriVals);
  end
  %******
  oo = floor(length(OriVals)/2);
  NPrefOri = [POri-oo-1,POri-oo,POri-oo+1];
  z = find(NPrefOri > length(OriVals));
  if ~isempty(z)
     NPrefOri(z) = NPrefOri(z) - length(OriVals);
  end
  z = find(NPrefOri <= 0);
  if ~isempty(z)
     NPrefOri(z) = NPrefOri(z) + length(OriVals);
  end
  %*******
  RTList = Info.SacOnList - Info.StimOnList;
  if (Info.MatchRT == 1) || ~isempty(Info.RTpass)
      BigList{1} = Info.MatchAttList; 
      BigList{2} = [ Info.MatchIgnAList ; Info.MatchIgnBList]; 
  else
      BigList{1} = Info.AttList; 
      BigList{2} = [ Info.IgnAList ; Info.IgnBList]; 
  end      
  
%   if isempty(RTpass)
%      BigList = cell(1,2);
%      BigList{1} = Info.AttList; 
%      BigList{2} = [ Info.IgnAList ; Info.IgnBList]; 
%   else
%      BigList = cell(1,2); 
%      attrts = RTList(Info.AttList);      
%      z = find( (attrts >= RTpass(1)) & (attrts < RTpass(2)) );
%      BigList{1} = Info.AttList(z);
%      uttrts = RTList(Info.IgnAList);
%      z1 = find( (uttrts >= RTpass(1)) & (uttrts < RTpass(2)) );
%      uttrts = RTList(Info.IgnBList);
%      z2 = find( (uttrts >= RTpass(1)) & (uttrts < RTpass(2)) );
%      BigList{2} = [ Info.IgnAList(z1); Info.IgnBList(z2) ];
%   end
%   %***********
%   %***********************************
%   %**** Resampling of trial distributions to match RTs
%   if (Info.MatchRT == 1)
%       OrigBigList = BigList;
%       BigList = cell(1,2);
%       RTList = Info.SacOnList - Info.StimOnList;
%       for n = 1:Info.MatchN
%         nBigList = spikestats.Match_Distribution_Lists(OrigBigList{1},OrigBigList{2},RTList,0.01);  % match in 10 ms bins, 0.01
%         BigList{1} = [BigList{1} ; nBigList{1}];
%         BigList{2} = [BigList{2} ; nBigList{2}];
%       end
%       %**** Note error bars will be artificially small, by sqrt(MatchN)
%   end
%   %*************
  
  Rast = cell(1,2);
  %*******
  tinfo.TimeLocks = TimeLocks;
  tinfo.AUC = cell(1,length(TimeLocks));
  tinfo.AUC_sem = cell(1,length(TimeLocks));
  tinfo.MI = cell(1,length(TimeLocks));
  tinfo.MIC = cell(3,length(TimeLocks));  % compute Nsmo = 0,1,2 (thus 3)
  tinfo.Fano = cell(1,length(TimeLocks));
  tinfo.RFano = cell(1,length(TimeLocks)); % rate-matched fano
  tinfo.Fstat = cell(1,length(TimeLocks));
  %****** for noise correlations on pairs later
  tinfo.SpkCnt = cell(2,length(TimeLocks));
  tinfo.SpkOri = cell(2,length(TimeLocks));
  %******
  for zk = 1:2  % zk is attention conditions
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
      %****** store raster *****
      staspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockWin(1); % in secs
      finspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockWin(2);
      z = find( (sp.st >= staspk) & (sp.st < finspk) );
      if ~isempty(z)
        difftime = (sp.st(z) - staspk );
        Rast{zk} = [Rast{zk} ; [(difftime+tinfo.LockWin(1)) (ones(size(difftime))*k)]];
      end
      %*************************
    end
    %*************
    tinfo.Event = Event;            % record what you time lock on
    tinfo.Rast{zk} = Rast{zk};
    tinfo.OriInd{zk} = OriDir{zk};  % trial list of orientation per trial (could do RT instead)
    tinfo.OriIndRT{zk} = OriRT{zk}; % need to figure this out
    VM = length(OriDir{zk});   % must give PSTH number of trials to fill rast
    [tt,uu,su] = spikestats.CompPSTH(tinfo.Rast{zk},tinfo.LockWin,VM,tinfo.Smooth);     
    tinfo.LockTT{zk} = tt;
    tinfo.LockUU{zk} = uu;
    tinfo.LockSU{zk} = su * sqrt(Info.MatchN);
    %******* compute AUC here
    tinfo.PrefOri = PrefOri;
    tinfo.NPrefOri = NPrefOri;
    [tt,auc,s_auc] = spikestats.CompCountSlider(tinfo.Rast{zk},PrefOri,NPrefOri,OriVals,OriInd{zk},tinfo.LockWin,tinfo.CountWin,1);
    tinfo.AucTT{zk} = tt;
    tinfo.AucUU{zk} = auc;
    tinfo.AucSU{zk} = s_auc * sqrt(Info.MatchN);
    %******* compute DSI over time
    DCountWin = 0.02;  % small counting bins for this, 20 ms width
    [tt,uu,su] = spikestats.CompCountSlider(tinfo.Rast{zk},[],[],OriVals,OriInd{zk},tinfo.LockWin,DCountWin,3);
    tinfo.DCountWin = DCountWin;
    tinfo.DSI_tt{zk} = tt;
    tinfo.DSI_uu{zk} = uu;
    tinfo.DSI_su{zk} = su * sqrt(Info.MatchN);  % maybe not expand here ...
    %**********
    for kk = 1:length(TimeLocks) 
       %******* then compute it for a specific interval defined by LockInt
       LockInt = TimeLocks{kk};
       OSpk = spikestats.CompCounts(LockInt(1),LockInt(2),OriInd{zk},tinfo.Rast{zk});
       [uval,sval,ci] = spikestats.Compute_AUC(PrefOri,NPrefOri,OriInd{zk},OSpk);
       tinfo.AUC{kk}(zk) = uval; 
       tinfo.AUC_sem{kk}(zk) = sval * sqrt(Info.MatchN);
       %*********
       cwin = LockInt(2)-LockInt(1);  % time window in seconds
       tinfo.MI{kk}(zk) = presac.Compute_mutualInformation(OriInd{zk},OSpk);
       %********* store counts for noise correlation computations later
       tinfo.SpkOri{zk,kk} = OriInd{zk};
       tinfo.SpkCnt{zk,kk} = OSpk;
       %***********
       for nn = 1:3
          NSmo = nn-1;
          tinfo.MIC{nn,kk}(zk) = presac.Compute_mutualInformation_circle(OriInd{zk},OSpk,NVals,NSmo,cwin);
       end
       ff = presac.Compute_fstat_fano(OriInd{zk},OSpk,cwin);
       tinfo.FanoSlope{kk}(zk) = ff(2);
       tinfo.Fstat{kk}(zk) = ff(1);
       tinfo.Fano{kk}(zk) = ff(3);
  %     tinfo.RFano{kk}(zk) = ff(4);  % rate-matched Fano factor 
                                     % accounts for non-linear var vs mean
    end
    %*********
  end % loop over zk conditions (att, ig)
  %****** compute a rate matched fano
  for kk = 1:length(TimeLocks)
       ff = presac.Compute_fano_rate_matched(tinfo.SpkOri{1,kk},tinfo.SpkCnt{1,kk},...
                                             tinfo.SpkOri{2,kk},tinfo.SpkCnt{2,kk},cwin);
       tinfo.RFano{kk}(1) = ff(1);
       tinfo.RFano{kk}(2) = ff(2);
       tinfo.RFanoSlope{kk}(1) = ff(3);
       tinfo.RFanoSlope{kk}(2) = ff(4);       
  end  
  %********
  if (0)
    disp('check on discrim stats');
    tinfo.AUC{1}
    tinfo.MI{1}
    tinfo.MIC{1,1}
    tinfo.Fano{1}
    tinfo.Fstat{1}
    input('done');
  end
  %**********************************
  %*********** save it to a new field *******
  Info = setfield(Info,FieldName,tinfo);
  %*********
  
return;

%******* below is the plot of stimulus locked PSTH
function plot_timelock_psth(Info,Field,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   subplot(H);
   tinfo = getfield(Info,Field);
   LockInt = tinfo.TimeLocks{1};
   spikeplot.plot_attention_traces(tinfo.LockTT,tinfo.LockUU,tinfo.LockSU,...
                                      tinfo.LockWin,LockInt,tinfo.Event);
   ylabel('Rate (sp/s)');
   title(sprintf('Name %s',Info.tagname));
   %****************
return;

%******* below is the plot of stimulus locked PSTH
function plot_timelock_auc(Info,Field,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   subplot(H);
   tinfo = getfield(Info,Field);
   LockInt = tinfo.TimeLocks{1};
   spikeplot.plot_attention_traces(tinfo.AucTT,tinfo.AucUU,tinfo.AucSU,...
                                      tinfo.LockWin,LockInt,tinfo.Event);
   ylabel('AUC');
   title(sprintf('AUC (%5.3f,%5.3f)',tinfo.AUC{1}(1),tinfo.AUC{1}(2)));
   %****************
return;

%******* below is the plot of stimulus locked PSTH
function plot_timelock_dsi(Info,Field,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   subplot(H);
   tinfo = getfield(Info,Field);
   LockInt = tinfo.TimeLocks{1};
   spikeplot.plot_attention_traces(tinfo.DSI_tt,tinfo.DSI_uu,tinfo.DSI_su,...
                                      tinfo.LockWin,LockInt,tinfo.Event);
   ylabel('DSI');
   %****************
return;
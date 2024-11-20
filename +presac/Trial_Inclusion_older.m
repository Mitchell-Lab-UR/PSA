function Info = Trial_Inclusion(Info,Exp,H)
%******* 
%******  function Info = Trial_Inclusion(Info,Exp,PlotIt)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          H - handle to plot window if provided, or if 1 create, else []
%***
%*** Outputs: Info - it updates fields in the Info struct to include 
%***                  information about trial inclusion parameters
%***

  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0) 
       hh = figure('position',[100 100 800 800]); 
       HA = subplot('position',[0.1 0.1 0.8 0.8]);
       plot_trial_inclusion(Info,HA);
       if (H == 2)  % store result to png and close
          z = getframe(hh);  % current figure
          uname = [Info.pathplot,filesep,'INC_',Info.tagname,'.png'];
          disp(sprintf('Storing image of Inclusion graph at %s',uname)); 
          imwrite(z.cdata,uname); % store image for later review
          close(hh);
       end
    else 
       if ~isempty(H)  % H should be a sub-plot handle otherwise
          plot_trial_inclusion(Info,H);
       end
    end
    return;
  end
  
  %****** parameters for determine trial correct to include in analyses
  TargRadius = 0.5;  % as a percentage of RF eccentricity, 
                     % what we countcas the saccade end point landing 
                     % close enough to the target, but should be relative 
                     % to the eccentricity (so a fraction, 0 to 1)
  AmpThresh = 0.5;  % fraction of saccade distance of total target ecc
                    % before we include that saccade far enough to target
  MinSacRT = 0.10;   % in seconds, min saccade onset from stim onset
  MaxAllowMicroSaccade = 0.5;  % largest micro-sacc allowed from stim onset to target saccade
  %*******************
  
  PlotDebug = 0;
  if (PlotDebug)
     HP = figure('position',[100,100,1200,400]);
  end
  
  %****************
  SacNum = 0;
  CorNum = 0;
  FixNum = 0;
  EarlyNum = 0;
  Tlist = [];  % list of trials with saccades correct
  TargList = [];  % list of which target he made a saccade towards
  TPosList = [];  % list of target choice, and x,y of sac endpoint
  TPosTrace = cell(1,1);
  PosCnt = 1;
  %***
  EPhysList = [];
  StimOnList = [];
  SacOnList = [];
  SacOffList = [];
  StimOri = [];
  %****
  SecondSac = 0;
  TrialCount = 0;
  for i = 1:size(Exp.D,1)
        if  strcmp(Exp.D{i}.PR.name(1:6),'FlagMo')
         
         bzz = [];       
         if isfield(Exp.D{i}.PR,'catchtrial')
             if (Exp.D{i}.PR.catchtrial == 1)  % do not include this trial
                 continue;
             end
         end
         
         if ~isfield(Exp.D{i}.PR,'target_null') || (Exp.D{i}.PR.target_null ~= 1)  % no blank target in RF location
               %********* 
               eyeSmo = Exp.D{i}.eyeSmo;
               eyeDat = Exp.D{i}.eyeData;
               slist = Exp.D{i}.slist;
               FixRadius = Exp.D{i}.P.fixWinRadius;  % within this bound for be fixating 
          if (Exp.D{i}.PR.error == 0) || (Exp.D{i}.PR.error == 4)  % only saccades to an aperture (no fix breaks, etc..) 
          
               TrialCount = TrialCount + 1;
               
               %***********
               stimori = Exp.D{i}.PR.targori;  % always in RF
               Tx = Exp.D{i}.PR.targ_x;
               Ty = Exp.D{i}.PR.targ_y;
               TargDist = norm([Tx(1),Ty(1)]);  % eccentricity of target
               BestTarg = NaN;
               BestEtx = NaN;
               BestEty = NaN;
               %********* for each trial, search for a saccade inside FixRadius at start
               %********* and in a target aperture at its end-point, and
               %********* was it a large saccade, not just stepping small
               if ~isempty(slist)
                  foundfix = 0;
                  foundsac = 0;
                  foundcor = 0;
                  for k = 1:size(slist,1)
                      stasac = slist(k,4); % integer start time of saccade
                      endsac = slist(k,5); % integer end time of saccade
                      stx = eyeSmo(stasac,2);
                      sty = eyeSmo(stasac,3);
                      etx = eyeSmo(endsac,2); 
                      ety = eyeSmo(endsac,3);
                      %******* did a saccade start in fix window?
                      fdist = norm([stx,sty]); 
                      if (foundfix == 0)
                          if (fdist < FixRadius)
                              foundfix = 1;
                          end
                      end
                      %****** if so, then did it end within an aperture
                      if (foundfix == 1)
                         edist = [];
                         for ik = 1:size(Tx,1)
                            edist(ik) = norm([Tx(ik),Ty(ik)]-[etx,ety]);
                         end
                         if (min(edist) < (TargDist * TargRadius) )  % first saccade into aperture
                            slist(k,7) = 1;  % mark the saccade
                            foundcor = 1;   % found a saccade to aperture
                            %****** compute fraction distance of saccade
                            sacamp = norm([(etx-stx),(ety-sty)]);
                            if (sacamp > (TargDist * AmpThresh) )
                               slist(k,7) = 2; % flag the saccade was large enough
                               foundsac = 1;
                               zz = find( edist == min(edist));
                               BestTarg = zz(1);   % which was target
                               BestTx = Tx(zz(1));
                               BestTy = Ty(zz(1));
                               BestEtx = etx;
                               BestEty = ety;
                               break;  % stop once the first is found
                            end
                            %**********
                         end
                      end
                  end  % end over loop to find big saccade, labeled slist(:,7) = 2
                  %******
                  FixNum = FixNum + foundfix;
                  CorNum = CorNum + foundcor;
                  SacNum = SacNum + foundsac;
                  %***********
                  
                  %**** that determines if the spatial aspect of the
                  %**** saccade was correct, but then two more question
                  %**** remain to include the trial:
                  %****   1) was the movement too early (relative to stim
                  %****                the stim onset)
                  %****   2) did a second saccade precede the big one
                  %**************
                  %*** Ques 1:  was it too early:
                  %*******
                  %*** Saccade Onset from Matlab states (state 5 in flight)
                  TooEarly = 0;  % start by assuming it is not too early
                  %**
                  MatStart = Exp.D{i}.eyeData(1,6);  
                  z = find( Exp.D{i}.eyeData(:,5) == 5);  % first moment Matlab detects leaving fix
                  if ~isempty(z)
                     msactime = Exp.D{i}.eyeData(z(1),1) - MatStart; 
                  else
                     % disp(sprintf('No saccade onset found trial %d',i)); 
                     msactime = NaN;
                  end
                  %**** Stim Onset from Matlab states
                  z = find( Exp.D{i}.eyeData(:,5) == 4);  %into state 4, stimuli appear
                  if ~isempty(z)
                    onstim = Exp.D{i}.eyeData(z(1),1) - MatStart;
                  else
                    onstim = NaN;
                  end
                  %*****
                  targsac = find( slist(:,7) == 2);
                  if ~isempty(targsac)
                      sac_ontime = slist(targsac(1),1);  % onset of that saccade
                      if ( (sac_ontime - onstim)  < MinSacRT)
                          TooEarly = 1;
                          EarlyNum = EarlyNum + 1;
                          continue;
                      end
                      sac_offtime = slist(targsac(1),2);  % offset of that saccade
                  else
                      sac_ontime = NaN;
                      sac_offtime = NaN;
                  end
                  
                  %********* then second question, are there any micro-saccades
                  %********* after stim onset and just before the main saccade
                  TooSecond = 0;
                  targsac = find( slist(:,7) == 2);
                  if ~isempty(targsac)
                      sac_ontime = slist(targsac(1),1);  % onset of that saccade
                      stasac = slist(targsac(1),4); % integer start time of saccade
                      %**** is there a micro-saccade just before
                      sac2nd = find( (slist(:,2) >= onstim) & (slist(:,2) < eyeSmo(stasac,1)) );
                      if (~isempty(sac2nd))
                           amp = norm( eyeSmo(slist(sac2nd,5),2:3) - eyeSmo(slist(sac2nd,4),2:3) );
                           if (amp > MaxAllowMicroSaccade)
                                TooSecond = 1;
                                SecondSac = SecondSac + 1;
                           end
                      end
                  end
                  
                  %******** final processing to list target
                  if (foundsac == 1) 
                       stasac = slist(targsac(1),4); % integer start time of saccade
                       endsac = slist(targsac(1),5); % integer start time of saccade
                       %*******
                       Tlist = [Tlist ; i];  % trial number to include
                       TargList = [TargList ; [BestTarg,BestTx,BestTy]]; % target choosen
                       %***** then target choose, and saccade end point
                       TPosList = [TPosList ; [BestTarg,BestEtx,BestEty]];
                       %******* store trace from before saccade
                       tzz = find( (eyeSmo(:,1) >= onstim) & ...
                                   (eyeSmo(:,1) < (eyeSmo(stasac,1)) ) );
                       bzz = find( (eyeSmo(:,1) >= onstim) & ...
                                       (eyeSmo(:,1) < (eyeSmo(endsac,1)) ) );
                       %*******
                       TPosTrace{PosCnt} = eyeSmo(tzz,1:3);
                       PosCnt = PosCnt + 1;
                       %********
                       if ~TooEarly && ~TooSecond && isfield(Exp.D{i},'START_EPHYS') 
                           EPhysList = [EPhysList ; 1];
                       else
                           EPhysList = [EPhysList ; 0];    
                       end
                       StimOnList = [StimOnList ; onstim];
                       SacOnList  = [SacOnList ; sac_ontime];
                       SacOffList = [SacOffList ; sac_offtime];
                       %****** stimulus ori in RF that trial
                       StimOri = [StimOri ; stimori];
                       %**************
                  end
                  %********   
               end
               
               %*************
               if (PlotDebug)
                   
                   figure(HP);
                   
                   Exp.D{i}.PR
                   [msactime,onstim,(msactime-onstim)]
                   [foundfix,foundcor,foundsac,TooEarly,TooSecond]
                                
                   
                   subplot('position',[0.1 0.1 0.3 0.8]); hold off;
                   plot(Tx,Ty,'bo'); hold on;
                   plot(0,0,'r+');
                   plot([-FixRadius,FixRadius],[0,0],'r-');
                   plot([0,0],[-FixRadius,FixRadius],'r-');
                   if ~isempty(bzz)
                      zz = find(eyeSmo(:,1) < eyeSmo(stasac,1));
                      plot(eyeSmo(zz,2),eyeSmo(zz,3),'g.-');             
                      zz = find(eyeSmo(:,1) > eyeSmo(endsac,1));
                      plot(eyeSmo(zz,2),eyeSmo(zz,3),'c.-');             
                      plot(eyeSmo(bzz,2),eyeSmo(bzz,3),'k.-');             
                   else
                      hh = plot(eyeSmo(:,2),eyeSmo(:,3),'k-');
                      set(hh,'Color',[0.7,0.7,0.7]);
                   end
                   axis([-10 10 -10 10]);
                   title(sprintf('SacRT %6.3f, Inc(%d,%d,%d)',...
                           (msactime-onstim),foundsac,TooEarly,TooSecond));
                   %******
                   subplot('position',[0.45 0.1 0.5 0.8]); hold off;
                   if ~isempty(bzz)
                      zz = find(( eyeSmo(:,1) > onstim) & (eyeSmo(:,1) < (eyeSmo(endsac,1)+0.1)) );
                      plot(eyeSmo(zz,1),eyeSmo(zz,2),'r-'); hold on;
                      plot(eyeSmo(zz,1),eyeSmo(zz,3),'b-');
                      axis([eyeSmo(zz(1),1) eyeSmo(zz(end),1) -10 10]);
                   else
                       plot(eyeSmo(:,1),eyeSmo(:,2),'r-'); hold on;
                       plot(eyeSmo(:,1),eyeSmo(:,3),'b-');    
                       axis([eyeSmo(1,1) eyeSmo(end,1) -10 10]);
                   end
                   for ik = 1:size(slist,1)
                       if (slist(ik,7) == 2) 
                          plot([slist(ik,1),slist(ik,1)],[-12,12],'m-');
                          plot([slist(ik,2),slist(ik,2)],[-12,12],'m-');
                       else
                          plot([slist(ik,1),slist(ik,1)],[-12,12],'k-');
                          plot([slist(ik,2),slist(ik,2)],[-12,12],'k-');                 
                       end
                   end
                   %******
                   input('view trial');
               end
               
               %********* store new saccade list
               Exp.D{i}.slist = slist;
               %*******************
          end
         end  % if trial was correct
     end
  end % if Flagmo trial
  %****** compute means of eye trace per target types
  tnum = length(unique(TPosList(:,1)));
  TPosTraceList = cell(1,tnum);
  for k = 1:size(TPosList,1)
     tg = TPosList(k,1);       
     if (EPhysList(k) == 1) 
       ux = nanmean(TPosTrace{k}(:,2));
       uy = nanmean(TPosTrace{k}(:,3));
       etx = TPosList(k,2);
       ety = TPosList(k,3);
       tx = TargList(k,2);
       ty = TargList(k,3);
       TPosTraceList{1,tg} = [TPosTraceList{1,tg} ; [ux,uy,etx,ety,tx,ty]];
     else
       TPosTraceList{1,tg} = [TPosTraceList{1,tg} ; [NaN,NaN,NaN,NaN,NaN,NaN]];
     end
  end
  %*********************
  
  %*************************
  disp('  ');
  disp(sprintf('Trial Count: %d ',TrialCount));
  disp(sprintf('Fixations %d, Correct %d', FixNum,CorNum));
  disp(sprintf('Saccades Flagged %d',SacNum));
  disp(sprintf('TooEarly %d SecondSac %d',EarlyNum,SecondSac));
  disp(sprintf('Final Included: %d  Ephys: %d',length(TargList),sum(EPhysList)));
  %*********************************************
  AttList = find( (TargList(:,1) == 1) & (EPhysList > 0) );
  IgnAList = find( (TargList(:,1) == 2) & (EPhysList > 0) );
  IgnBList = find( (TargList(:,1) == 3) & (EPhysList > 0) );
  if (length(IgnBList) > length(IgnAList) )
      tmp = IgnBList;
      IgnBList = IgnAList;
      IgnAList = tmp;
      %*** want more here to classify same or opposite hemifield?  Not now
  end
  
  disp(sprintf('Presaccadic conds (%d,%d,%d)',length(AttList),length(IgnAList),length(IgnBList)));
  %**********
  Info.FixNum = FixNum;
  Info.CorNum = CorNum;
  Info.SacNum = SacNum;
  Info.EarlyNum = EarlyNum;
  Info.SecondSac = SecondSac;
  Info.Tlist = Tlist;  % list of trials with saccades correct
  Info.TargList = TargList;  % list of which target he made a saccade towards
  Info.TPosList = TPosList;  % list of target choice, and x,y of sac endpoint
  Info.TPosTrace = TPosTrace; % traces of eye position prior to saccades
  Info.TPosTraceList = TPosTraceList;  % list of mean eye pos bef saccade
  Info.EPhysList = EPhysList;  % list if trial has Ephys to use
  Info.StimOnList = StimOnList; % Stimulus onset times
  Info.SacOnList = SacOnList;   % Saccade onset times
  Info.SacOffList = SacOffList; % Saccade offset times
  Info.AttList = AttList;   % Ephys trials with target in RF
  Info.IgnAList = IgnAList; % Ephys trials target outside RF, most trials
  Info.IgnBList = IgnBList; % Ephys trials target outside RF, fewer trials
  Info.StimOri = StimOri;   % Direction of stimulus inside RF  
  
return;

%******* below is the plot function in the same file
function plot_trial_inclusion(Info,H)
   %******* This plots 1) saccade end points around targets
   %***                2) the eye traces from stim onset to saccade onset
   %***                     per color by target
   %***                3) the mean and sem of pre-sac trace position
   
   subplot(H);
   showrawtrace = 1;  % show pre-sac traces
   colo = 'rbcg';
   for k = 1:size(Info.TPosList,1)  
     tg = Info.TPosList(k,1);
     if (showrawtrace) && (Info.EPhysList(k)==1)
        h = plot(Info.TPosList(k,2),Info.TPosList(k,3),[colo(tg),'o']); hold on;
        set(h,'Markersize',4);
        %*** plot eye position before the saccade
        h = plot(Info.TPosTrace{k}(:,2),Info.TPosTrace{k}(:,3),[colo(tg),':']); hold on;
     end
     %******
  end
  %****** get mean eye position before saccade 
  for k = 1:length(Info.TPosTraceList)
      ux = nanmean(Info.TPosTraceList{1,k}(:,1));
      uy = nanmean(Info.TPosTraceList{1,k}(:,2));
      h = plot(ux,uy,[colo(k),'o']); 
      set(h,'Markersize',10); hold on;
      set(h,'Linewidth',3);
      tux = nanmean(Info.TPosTraceList{1,k}(:,3));
      tuy = nanmean(Info.TPosTraceList{1,k}(:,4));
      h = plot(tux,tuy,[colo(k),'s']); 
      set(h,'Markersize',10);
      set(h,'Linewidth',3);
      tx = nanmean(Info.TPosTraceList{1,k}(:,5));
      ty = nanmean(Info.TPosTraceList{1,k}(:,6));
      h = plot(tx,ty,'k+'); 
      set(h,'Markersize',10);
      set(h,'Linewidth',3);
  end
  %**************
  h = plot(0,0,'k+');
  set(h,'Markersize',10);
  maxo = max(max(abs(Info.TPosList(:,2:3))));
  axis([-maxo maxo -maxo maxo]);

return;
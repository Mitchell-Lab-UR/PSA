function Info = Compute_Flash_Response(Info,Exp,Unit,H,FieldName,Event,WinLock,UseOri,NoiseType)
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
%***          Event - if 1, flash onset 
%***                  if 2, flash offset
%***          WinLock - window of analyses
%***          UseOri - sort raster based on orientation
%***          NoiseType - type of flash
%***
%*** Outputs: Info - it computes the flash-response PSTH and rasters

  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh = figure('position',[200 70 800 800]); 
       HA = subplot('position',[0.075 0.60 0.375 0.35]);  % PSTH
    else
        if isempty(H)
          return;
        else
          if (length(H) ~= 2)
              disp('Error with subplot pass to Compute_DSI, requires two subplots');
              return;
          else
              HA = H{1};
          end
        end
    end
    %*******
    plot_timelock_psth(Info,FieldName,HA);
    %********
    if isfloat(H)  && (H == 2)
       z = getframe(hh);  % current figure
       uname = [Info.pathplot,filesep,FieldName,'_FlashResponse_',Info.tagname,'.png'];
       disp(sprintf('Storing image of TimeLock graph at %s',uname)); 
       imwrite(z.cdata,uname); % store image for later review
       close(hh);
    end   
    return;
  end
  
  %****** parameters for computing stimlock response 
  tinfo = []; % new struct for this local 'FieldName'
  tinfo.LockWin = WinLock;    % plot stim locked firing rate PSTH 
  tinfo.LockInt = [0.05 0.10];
  if isempty(UseOri)
      tinfo.Smooth = 2;
      UseOri = 0;
  else
      tinfo.Smooth = 5;           % Gaussian sig for smooth PSTH plots
  end
  tinfo.NoiseType = NoiseType;  % 3 or 6
  %*******************
  
  %****** compute the DSI ******
  %**** make an Uber Unit with all array spikes
  if (1)
    sp.st = [];
    for k = 1:length(Exp.sp)
      sp.st = [sp.st ; Exp.sp{k}.st];
    end
  else
    sp = Exp.sp{Unit};
  end
  %*******
 
  [FlashHits,FlashTlist,FlashOri] = getFlashEvents(Exp,1,NoiseType);  
  [FlashHits2,FlashTlist2,FlashOri2] = getCSDEvents(Exp,NoiseType);
  
  if (0)  % looks like codes come back as identical
    figure(10);
    subplot(2,2,1);
    plot(FlashHits,FlashHits2,'ko');
    subplot(2,2,2);
    plot(FlashTlist,FlashTlist2,'ko');
    input('hold');
  end
  if (1)
      FlashHits = FlashHits2;
      FlashTlist = FlashTlist2;
      FlashOri = FlashOri2;
  end
  
  tinfo.FlashHits = FlashHits;
  tinfo.FlashTlist = FlashTlist;
  tinfo.FlashOri = FlashOri;
  if ~isempty(FlashHits)   
      %************************
      Rast = [];
      Ori = []; % integer orientation values per trial
      %***** for reference, what are the set of orientations
      disp('Collecting spike counts per orientation for stimulus tuning ...');
      for k = 1:length(FlashTlist)
          trial = FlashTlist(k); % reference to Exp trial
          tr = k;  % reference in trial counts directly
          %****** get stim spike count on trial
          tlock = FlashHits(tr); 
          %****** store raster *****
          staspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockWin(1); % in secs
          finspk = Exp.D{trial}.START_EPHYS + tlock + tinfo.LockWin(2);
          z = find( (sp.st >= staspk) & (sp.st < finspk) );
          if ~isempty(z)
            difftime = (sp.st(z) - staspk );
            Rast = [Rast ; [(difftime+tinfo.LockWin(1)) (ones(size(difftime))*k)]];
            Ori = [Ori ; FlashOri(tr)];
          end
          %*************************
      end 
      tinfo.Rast = Rast;
      tinfo.Ori = Ori;
      tinfo.Event = Event;            % record what you time lock on
      tinfo.UseOri = UseOri;          % use orientation to raster index
      VM = length(Ori);   % must give PSTH number of trials to fill rast
      if ~isempty(Rast)
        [tt,uu,su] = spikestats.CompPSTH(Rast,tinfo.LockWin,VM,2); % tinfo.Smooth);     
        tinfo.LockTT = tt;
        tinfo.LockUU = uu;
        tinfo.LockSU = su;
      else
        tinfo.LockTT = [];
        tinfo.LockUU = [];
        tinfo.LockSU = [];         
      end
      %*************
  end
  
  %*********** save it to a new field *******
  Info = setfield(Info,FieldName,tinfo);
  %*********
  
return;

%******* below is the plot of stimulus locked PSTH
function plot_timelock_psth(Info,Field,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   subplot(H);
   tinfo = getfield(Info,Field);
   if isfield(tinfo,'LockTT')
       if ~isempty(tinfo.LockTT)
           plot(tinfo.LockTT,tinfo.LockUU,'k-'); hold on;
           plot(tinfo.LockTT,tinfo.LockUU+(2*tinfo.LockSU),'k-'); hold on;
           plot(tinfo.LockTT,tinfo.LockUU-(2*tinfo.LockSU),'k-'); hold on;
       end
   end
   ylabel('Rate (sp/s)','FontSize',15');
   title(sprintf('Name %s',Info.tagname),'FontSize',18);
   %*************
return;

function [FlashHits,FlashTlist,FlashOri] = getFlashEvents(Exp,Event,NoiseType)
  
  %****** search for trial events in the Exp file
  FlashHits = [];
  FlashOri = [];  % set as 1 if only a contrast flash
  FlashTlist = []; % list of trials
  %***********
  for i = 1:size(Exp.D,1)
     if  strcmp(Exp.D{i}.PR.name(1:6),'Forage')
        if (Exp.D{i}.PR.noisetype == NoiseType)  % CSD flashes
            %*******************
            MatStart = Exp.D{i}.eyeData(6,1);  
            nhist = Exp.D{i}.PR.NoiseHistory;
            for k = 2:length(nhist)
                if (Event == 1)  % flash onsets
                  if (nhist((k-1),2)==0) && (nhist(k,2)>0)  % flash onset
                     ontime = nhist(k,1) - MatStart;  % time relative trial start
                     FlashTlist = [FlashTlist ; i];  % trial start ref
                     FlashHits = [FlashHits ; ontime]; % relative start ref
                     FlashOri = [FlashOri ; nhist(k,2)];  % ori index
                  end
                end
                %*****
                if (Event == 2)  % flash offsets
                  if (nhist((k-1),2)>0) && (nhist(k,2)==0)  % flash offset
                     ontime = nhist(k,1) - MatStart;  % time relative trial start
                     FlashTlist = [FlashTlist ; i];  % trial start ref
                     FlashHits = [FlashHits ; ontime]; % relative start ref
                     FlashOri = [FlashOri ; nhist(k,2)];  % ori index
                  end
                end
            end
            %*******************
        end
     end
  end
 return;
  
function [FlashHits,FlashTlist,FlashOri] = getCSDEvents(Exp,which_noisetype)
   
%function eventTimes = getCSDEventTimes(Exp, which_noisetype)
% gets the event times of the current source density trials
% Inputs:
%   Exp              [struct] - Exp struct from io.dataFactoryGratingSubspace
% 
% jfm wrote it 2020
% ghs edit it 2020

if nargin < 2
  which_noisetype = 3;
end

sanitycheck = 0; % check spike histogram to verify CSD activity

Trials = length(Exp.D);
CSD = [];
CSD.Trials = [];
for k = 1:Trials
   ProtoName = Exp.D{k}.PR.name;
   %ProtoName
   type = 0;
   if (strcmp(ProtoName,'ForageProceduralNoise'))  % for Jake
       type = 1;
   end
   if (strcmp(ProtoName,'Forage'))   % for Shanna
       type = 2;
   end
  % type
   if (type > 0)
       NoiseType = Exp.D{k}.PR.noisetype;
       %disp(sprintf('Trial(%d); %s  NoiseType(%d)',k,ProtoName,NoiseType));
       if (type == 1) && (NoiseType == 3)
           CSD.Trials = [CSD.Trials ; [k 1]];  % 1 indicate contrast onset
       end
       if (type == 2) && ( (NoiseType == 3) || (NoiseType == 6) ) 
           if (NoiseType == 3) && (which_noisetype==3)
              CSD.Trials = [CSD.Trials ; [k 2]];
           end
           if (NoiseType == 6) && (which_noisetype==6)
              CSD.Trials = [CSD.Trials ; [k 3]];
           end
       end
   end
end

FlashTlist = []; 
FlashHits = [];
FlashOri = [];

CSD.Onsets = [];
CSD.Types = [];
CSD.MoDir = [];
CSD.Offsets = [];
CSD.Onsets_Ephys = [];
CSD.Offsets_Ephys = [];
NTrials = length(CSD.Trials);
for k = 1:NTrials
    kk = CSD.Trials(k,1);
    type = CSD.Trials(k,2);
    NoHist = Exp.D{kk}.PR.NoiseHistory;
    MatStart = Exp.D{kk}.eyeData(6,1);  % Matlab start
    
    %*** Noise History is time in column 1, and contrast (0 or 1) in col 2
    %** find all the Onsets, as transition from 0 to 1 in column 2
    for i = 2:size(NoHist,1)
       if (NoHist(i-1,2) == 0) && (NoHist(i,2) >= 1)   % 
           CSD.Onsets = [CSD.Onsets ; NoHist(i,1)];  % store Mat time
           CSD.Types = [CSD.Types ; type];  % 1 - contrast (Jake), 2 - contrast (Shanna), 3 - motion (Shanna)
           CSD.MoDir = [CSD.MoDir ; NoHist(i,2)];
           
           %*** DOUBLE TEETH PROBLEM MAY ARISE FROM USING
           %*** 'STARTCLOCKTIME' ... be safe, use Matlab Trial Start
           %*****
           %******* convert to Ephys time per trial start-clocks
           % if isfield(Exp.D{kk},'STARTCLOCKTIME')
           %    tt = tt - Exp.D{kk}.STARTCLOCKTIME;  % Jake - 0 from start of trial
           % else
           %    tt = tt - Exp.D{kk}.eyeData(6,1);    % Shanna - start of trial in mat time
           % end
           
           tt = NoHist(i,1) - MatStart + Exp.D{kk}.START_EPHYS;     % time from start of trial in ephys
           CSD.Onsets_Ephys = [CSD.Onsets_Ephys ; tt];
           
           %*************
           FlashTlist = [FlashTlist ; kk];
           FlashHits = [FlashHits ; (NoHist(i,1)-MatStart)];
           FlashOri = [FlashOri ; NoHist(i,2)];
           %************
           
           % search for corresponding offset, if trial ends, insert a NaN
           for j = (i+1):size(NoHist,1)
             testOff = 0;
             if (type == 3)
                 testOff = (NoHist(j-1,2) >= 0) && (NoHist(j,2) < 0);
             else
                 testOff = (NoHist(j-1,2) >= 1) && (NoHist(j,2) == 0);
             end
             if (testOff)
               CSD.Offsets = [CSD.Offsets ; NoHist(j,1)];  % store Mat time
               
               %*** DOUBLE TEETH PROBLEM MAY ARISE FROM USING
               %*** 'STARTCLOCKTIME' ... be safe, use Matlab Trial Start
               %*****
               %******* convert to Ephys time per trial start-clocks
               % if isfield(Exp.D{kk},'STARTCLOCKTIME')
               %  tt = tt - Exp.D{kk}.STARTCLOCKTIME;  % 0 from start of trial
               % else
               %  tt = tt - Exp.D{kk}.eyeData(6,1);    % Shanna - start of trial in mat time 
               % end
               
               tt = NoHist(j,1) - MatStart + Exp.D{kk}.START_EPHYS;     % time from start of trial in ephys
               CSD.Offsets_Ephys = [CSD.Offsets_Ephys ; tt];
               break;
             end
           end
           if (j >= size(NoHist,1))  % Offset occured after trial end
               CSD.Offsets = [CSD.Offsets ; NaN];
               CSD.Offsets_Ephys = [CSD.Offsets_Ephys ; NaN];
           end
       end 
    end
    %*** sanity check result looks right per trial
    if (sanitycheck == 1)
        figure(10); hold off;
        plot(NoHist(:,1),NoHist(:,2),'k-'); hold on;
        % plot(CSD.Onsets,(0.5*ones(size(CSD.Onsets))),'rx');
        xlabel('Time (secs)');
        %input('check');
    end
end   
return;


  


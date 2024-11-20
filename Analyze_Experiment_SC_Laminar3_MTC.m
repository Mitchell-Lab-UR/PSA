function Analyze_Experiment_SC_Laminar3_MTC(FileTag,ProcPath,StorePath,PlotPath,PlotIt,UseFoveal)
%***
%*** function Analyze_Experiment(FileTag,ProcPath,StorePath,PlotPath,PlotIt)
%***
%*** inputs: FileTag - tag of session to analyze (Exp file in Proc dir)
%***         ProcPath - folder with processed Exp files
%***         StorePath - folder where to store results of analyses
%***         PlotPath - folder where to plot results stored as .png
%***         PlotIt - show results as you go
%***
%*** general organization:
%***     analysis will create many partial results to be stored
%***     store them as 'FileTag_U%d_name.mat' where name is the analysis type
%***     
%***     each analysis produces its own output product, which is stored
%***     that product can be loaded and shown in a plot
%***     so for each analysis, define 1) data struct returned
%***                                  2) analysis function and 
%***                                  3) plot function

%%
%*** designed to take inputs from a batch file, but has defaults if none
RunByParts = 0;
PlotIt = 0;
PlotIt2 = 0;
% Animal='Processed_Sprout';
% User='abucklaew';

if RunByParts || isempty(FileTag)
   %FileTag = 'Sprout_251022'; 
end

if RunByParts || isempty(ProcPath)
    PROCESSED_DATA_DIR = 'C:\Users\abucklaew\Box\MTC\Processed_Sprout';
   % PROCESSED_DATA_DIR = 'C:\Users\amybu\Box\MTC\Processed_Sprout';
   %PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Processed_Sprout';
else
   PROCESSED_DATA_DIR = ProcPath;
end

if RunByParts || isempty(StorePath)
  % STORE_DATA_DIR = 'C:\Users\amybu\Box\MTC\Info_Sprout'; 
   STORE_DATA_DIR = 'C:\Users\abucklaew\Box\MTC\Info_Sprout_MT'; 
  %STORE_DATA_DIR = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Info_Sprout';
else
   STORE_DATA_DIR = StorePath; 
end

if RunByParts || isempty(PlotPath)   
    %PLOT_DATA_DIR = 'C:\Users\amybu\Box\MTC\Plot_Info';
    PLOT_DATA_DIR = 'C:\Users\abucklaew\Box\MTC\Plot_Info';
   %PLOT_DATA_DIR = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Plot_Info';
else
   PLOT_DATA_DIR = PlotPath; 
end
%*************

%% Grab Exp File that was previously imported
ExpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'.mat'];
disp(sprintf('Loading processed Exp file, %s',FileTag));
load(ExpFile);  % struct Exp will be loaded
disp('File loaded');

%% ******* Load in Waveform information (added by Jude 7/11)  NEW CODE (7/11)
if isfield(Exp,'osp')
  [SpikeDurations,SpikeWaves,SpikeChan] = OspToWaves(Exp.osp);
else
  SpikeDurations = NaN;
end
%******************

%******* Behavior can be computed outside the unit loop and repeated
if (1)    
    Info.tagname = FileTag;
    %******* Enforce matching of RT distributions
    if (1)
       Info.MatchRT = 1;
       Info.MatchN = 8;
       % NOTE: Stim transient is typically 40-80ms, and presac -50 to 25
       %       So setting min RT at 130 insures no stim transient in presac
       %       Which is cleaner, plus may reduce variance (though fewer
       %       trials will be included, lose about 20-30% possibly
       Info.RTpass = [0.120 Inf];  % only include RTs after 130 ms lists
    else
       Info.MatchRT = 0;
       Info.MatchN = 1;    
       Info.RTpass = [];
    end
    
    %***** compute trial inclusion list from criteria, pass it Exp and Info
%     if (UseFoveal == 1)
%        if (fname(1) == 'E')
%           FOVANG = 180;
%           FOVERR = 50;
%        else
%           FOVANG = 120;
%           FOVERR = 50; 
%        end
%        Info = presac.Trial_Inclusion_SC_FOV(Info,Exp,0,FOVANG,FOVERR);
%        if isempty(Info)
%            continue;
%        end
%        if (PlotIt)
%           presac.Trial_Inclusion_SC_FOV(Info,[],PlotIt,FOVANG,FOVERR);
%        end       
%     else
    
    if (1)
       % set FOVANG = [] if not using that option
       Info = presac.Trial_Inclusion_SC_FOV(Info,Exp,0,[],[]);
       if (PlotIt)
          presac.Trial_Inclusion_SC_FOV(Info,[],PlotIt,[],[]);
       end
    end
    
    %*********LEAVE OUT PFR COMPUTATIONS FOR NOW *************
    % plotsteps = 0;
    % trialtype = NaN;
    % Info = presac.Compute_PFR_Model2d(Info,Exp,0,plotsteps,20);
    % if (PlotIt)
    %       presac.Compute_PFR_Model2d(Info,[],PlotIt,0,20);
    % end
    
    %***** compute basic saccade stats ... yes need for foraging figure
    Info = presac.Compute_SaccadeFix(Info,Exp,0);
    if (PlotIt)
         presac.Compute_SaccadeFix(Info,[],PlotIt);
    end
    
    KeepInfo = Info;  % duplicate for units coming
end


%******************
%% Compute trial inclusion and spiking stats per unit
for UNIT = 1:max(length(Exp.sp),1)
%UNIT
% for UNIT = MyUnits %****
    %*********
    fname = sprintf('%s_U%d',FileTag,UNIT);
    unitname = [STORE_DATA_DIR,filesep,fname,'.mat']; 
    
    %% ******* If you wish to expedite a routine and append to Info
    if (0)     % load pre-saved file to expedite processing
       if exist(unitname)
          disp(sprintf('Loading previous Info %s to append it',unitname));
          load(unitname);
       else
          continue;
       end

       %****** compute PSTH per condition, Stim Onset locked tuning
%        if (Info.Included == 1)
%           Event = 1; % stim onset
%           WinLock = [-0.1,0.2]; 
%           TimeLock = [0.04,0.09];
%           UseOri = []; % if UseOri = [], then will make smoothing 2ms instead of 5ms
%           %*******
%           Info = presac.Compute_TimeLock(Info,Exp,UNIT,0,'StimOn2',...
%                                           Event,WinLock,TimeLock,UseOri); 
%           if (PlotIt)
%                  presac.Compute_TimeLock(Info,[],UNIT,PlotIt,'StimOn2',...
%                                       Event,WinLock,TimeLock,UseOri);  %if PlotIt=2, writes .png output    
%           end
%           %**********
%           Event = 1; % flash onset
%           WinLock = [-0.1,0.2]; 
%           UseOri = []; % use orientation to index raster
%           %*****
%           NoiseType = 3;
%           Info = presac.Compute_Flash_Response(Info,Exp,UNIT,0,'FlashResp32',Event,WinLock,UseOri,NoiseType); 
%           if (PlotIt)
%               presac.Compute_Flash_Response(Info,[],UNIT,PlotIt,'FlashResp32',Event,WinLock,UseOri,NoiseType);  %last input 2, writes .png output    
%           end
%           %*******
%           NoiseType = 6;
%           Info = presac.Compute_Flash_Response(Info,Exp,UNIT,0,'FlashResp62',Event,WinLock,UseOri,NoiseType); 
%           if (PlotIt)
%              presac.Compute_Flash_Response(Info,[],UNIT,PlotIt,'FlashResp62',Event,WinLock,UseOri,NoiseType);  %last input 2, writes .png output    
%           end
          
          
           
%           Info.SacImage = Forage.Compute_Saccade_Lock(Info,Exp,UNIT,0,1); 
%           if (PlotIt)
%             if ~isempty(Info.SacImage)
%                 Forage.Compute_Saccade_Lock(Info.SacImage,[],UNIT,PlotIt,1);  %last input 2, writes .png output  
%                 input('view back image modulation');
%                 close all;
%             end
%           end
          
          %*** STORE results to info struct, one info per neuron in database
          save(unitname,'Info');
          disp('saved file')
          if (PlotIt == 2) || (PlotIt2 > 0)
            input('check unit');
            close all;
          end
          %*** continue
          clear Info;
       end
       %continue  % skip other analyses and just go to next one
       
    %end
    %************
    
    %%
    Info = KeepInfo;  % copy over primary element with behavior and trial info
    Info.tagname = fname;
    Info.pathname = STORE_DATA_DIR;
    Info.pathplot = PLOT_DATA_DIR;
    if isfield(Exp,'sp') && ~isempty(Exp.sp) && isfield(Exp.sp{UNIT},'iso')
       if ~isfield(Exp.sp{UNIT},'nulldist')
           Info.isolation = [Exp.sp{UNIT}.iso NaN NaN];
           disp(sprintf('Note for %s, it is missing spike Nulldist',fname));
           %input('Press return to continue anyway');
       else
           Info.isolation = [Exp.sp{UNIT}.iso,...
                  Exp.sp{UNIT}.nulldist,Exp.sp{UNIT}.neardist];
       end
    else
       Info.isolation = NaN;
    end 
    
    
    if isfield(Exp,'PHYClu') && isfield(Exp,'osp')
       Info.depth = (Exp.osp.clusterDepths(1,UNIT));
       if (size(Exp.PHYClu,1) < UNIT)
           disp(sprintf('UNABLE TO RECOVER ISOLATION, %s',fname));
           Info.isolation = [0 NaN NaN];  % default noise
       else
           if strcmp(Exp.PHYClu.group(UNIT),'good')  %% what to set this as? 
               Info.isolation = [4 NaN NaN];
           end
           if strcmp(Exp.PHYClu.group(UNIT),'mua')  %% what to set this as? 
               Info.isolation = [1 NaN NaN];
           end
           if strcmp(Exp.PHYClu.group(UNIT),'noise')  %% what to set this as? 
               Info.isolation = [0 NaN NaN];
           end
       end
    else 
       if isfield(Exp,'osp')
          Info.depth = (Exp.osp.clusterDepths(1,UNIT));
       else
          Info.depth = NaN;
       end
    end
    
    %********** added by Jude, 7/11/2021
    if isfield(Exp,'osp')
      Info.duration = SpikeDurations(UNIT);
      Info.waveform = SpikeWaves(UNIT,:);  % bug was here
      Info.channel = SpikeChan(UNIT); % added by Shanna 7/20/21 to get Shank info
      if SpikeChan(UNIT) > 32
        Info.shank = 2;
      else 
        Info.shank = 1; 
      end
    else
      Info.waveform = [];
      Info.duration = NaN;   % possible to extract actually
      Info.channel = NaN;
      Info.shank = 0;
    end
    %**************
    
    %% ****
%     for UNIT = 1:max(length(Exp.sp),1)
%     UNIT
%     PlotIt = 1;
    
    % Compute saccade modulation for BackImage trials, Noisetype = Nan (not use)
    UseBackImage = 1;
    Info.SacImage = Forage.Compute_Saccade_Lock(Info,Exp,UNIT,0,UseBackImage); 
    if (PlotIt)
        if ~isempty(Info.SacImage)
          Forage.Compute_Saccade_Lock(Info.SacImage,[],UNIT,PlotIt,UseBackImage);  %last input 2, writes .png output  
          input('view back image modulation');
          close all;
        end
    end
    
    %end
    
    %% ********
    %****** compute DSI metric forunit
    Info = presac.Compute_DSI(Info,Exp,UNIT,0); 
    if (PlotIt2)
        presac.Compute_DSI(Info,[],UNIT,PlotIt2);  %last input 2, writes .png output    
    end
    
    %% ****** compute Hartley metric for each unit
    %*****Added by Amy on 11/22/21
    if (PlotIt2)==0;
        Info = presac.Compute_Hartley2(Info,Exp,UNIT,0); 
    else
        Info = presac.Compute_Hartley2(Info,Exp,UNIT,PlotIt2);     
    end
    %input('check');
   
    %%
    %****** compute Flash Response to contrast onset (if CSD was done, laminar)
    if (Info.Included == 1)  % requires neuron had visual response at least
      Event = 1; % flash onset
      WinLock = [-0.0,0.8]; 
      UseOri = 1; % use orientation to index raster
      %*****
      NoiseType = 3;
      Info = presac.Compute_Flash_Response(Info,Exp,UNIT,0,'FlashResp3',Event,WinLock,UseOri,NoiseType); 
      if (PlotIt)
        presac.Compute_Flash_Response(Info,[],UNIT,PlotIt,'FlashResp3',Event,WinLock,UseOri,NoiseType);  %last input 2, writes .png output    
      end
      %*******
      % NoiseType = 6;
      %Info = presac.Compute_Flash_Response(Info,Exp,UNIT,0,'FlashResp6',Event,WinLock,UseOri,NoiseType); 
      %if (PlotIt)
      %  presac.Compute_Flash_Response(Info,[],UNIT,PlotIt,'FlashResp6',Event,WinLock,UseOri,NoiseType);  %last input 2, writes .png output    
      %end
    end
    
    %% 
    %****** compute behavioral performance
    %****** compute PSTH per condition, Stim Onset locked tuning
    if (Info.Included == 1)
       Event = 1; % stim onset
       WinLock = [-0.1,0.2]; 
       TimeLock = [0.04,0.09];
       UseOri = 1; % use orientation to index raster
       %*******
       Info = presac.Compute_TimeLock(Info,Exp,UNIT,0,'StimOn',...
                                          Event,WinLock,TimeLock,UseOri); 
       if (PlotIt)
              presac.Compute_TimeLock(Info,[],UNIT,PlotIt,'StimOn',...
                                      Event,WinLock,TimeLock,UseOri);  %if PlotIt=2, writes .png output    
       end           
    end
    
    %****** compute PSTH per condition, Saccade Onset locked tuning
    if (Info.Included == 1)
       Event = 2; % sac onset
       WinLock = [-0.20,0.05];
       TimeLock = [-0.03,0.03];
       UseOri = 0; % use orientation to index raster
       %******
       Info = presac.Compute_TimeLock(Info,Exp,UNIT,0,'SacOn',...
                                          Event,WinLock,TimeLock,UseOri); 
       if (PlotIt2)
                   presac.Compute_TimeLock(Info,[],UNIT,PlotIt2,'SacOn',...
                                      Event,WinLock,TimeLock,UseOri);  %if PlotIt=2, writes .png output        
       end           
    end
         
    %****** compute AUC and Mutal Info, stimulus onset *****
    if (Info.Included == 1)
       Event = 1; % stim onset
       WinLock = [-0.05,0.20];
       TimeLocks = cell(1,1);
       TimeLocks{1} = [0.04,0.09];
       Info = presac.Compute_Discrim(Info,Exp,UNIT,0,'AUC_StimOn',...
                                          Event,WinLock,TimeLocks); 
       if (PlotIt)
              presac.Compute_Discrim(Info,[],UNIT,PlotIt,'AUC_StimOn',...
                                          Event,WinLock,TimeLocks);  %if PlotIt=2, writes .png output    
       end         
       %********  
    end

    %****** compute AUC, presaccadic period *****
    if (Info.Included == 1)
       Event = 2; % sac onset
       WinLock = [-0.20,0.10];
       TimeLocks = cell(1,1);
       TimeLocks{1} = [-0.03,0.03];
       Info = presac.Compute_Discrim(Info,Exp,UNIT,0,'AUC_SacOn',...
                                          Event,WinLock,TimeLocks); 
       if (PlotIt)
             presac.Compute_Discrim(Info,[],UNIT,PlotIt2,'AUC_SacOn',...
                                          Event,WinLock,TimeLocks);  %if PlotIt=2, writes .png output    
       end           
    end  

    %****** compute PSTH per condition, Saccade Onset locked tuning
    if (Info.Included == 1)
       Event = 3; % sac offset
       WinLock = [-0.10,0.15];
       TimeLock = [-0.06,0.04];  %post-saccadic before new visual input
       UseOri = 0; % use orientation to index raster
       %******
       Info = presac.Compute_TimeLock(Info,Exp,UNIT,0,'SacOff',...
                                          Event,WinLock,TimeLock,UseOri); 
       if (PlotIt2)
                   presac.Compute_TimeLock(Info,[],UNIT,PlotIt2,'SacOff',...
                                      Event,WinLock,TimeLock,UseOri);  %if PlotIt=2, writes .png output        
       end           
    end
    %*****************

    if (Info.Included == 1)
       Event = 3; % sac offset
       WinLock = [-0.10,0.15];
       TimeLocks = cell(1,1);
       TimeLocks{1} = [-0.06,0.04];
       Info = presac.Compute_Discrim(Info,Exp,UNIT,0,'AUC_SacOff',...
                                          Event,WinLock,TimeLocks); 
       if (PlotIt)
             presac.Compute_Discrim(Info,[],UNIT,PlotIt2,'AUC_SacOff',...
                                          Event,WinLock,TimeLocks);  %if PlotIt=2, writes .png output    
       end           
    end  
    
    %*** STORE results to info struct, one info per neuron in database
    save(unitname,'Info');
    disp('saved file')
    if (PlotIt == 2) || (PlotIt2 > 0)
        input('check unit');
        close all;
    end
    %*** continue
    clear Info;
    %*******************
    
end


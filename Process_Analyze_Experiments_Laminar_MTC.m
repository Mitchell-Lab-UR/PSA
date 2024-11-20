%%  ********* This is a script to run analysis on multiple experiments
%***
%***
if (0) % Jude's Laptop
   ProcPath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Processed_MTC';
   StorePath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Info_MTC';
   PlotPath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Plot_Info';
else   % Amy's CPU
%     ProcPath = 'C:\Users\amybu\Box\MTC\Processed_MTC';
%     StorePath = 'C:\Users\amybu\Box\MTC\Info_MTC';
%     PlotPath = 'C:\Users\amybu\Box\MTC\Plot_Info';
    ProcPath = 'C:\Users\abucklaew\Box\MTC\Processed_MTC';
    StorePath = 'C:\Users\abucklaew\Box\MTC\Info_MTC';
    PlotPath = 'C:\Users\abucklaew\Box\MTC\Plot_Info';
    
end
%******** Folder list for updates
FolderList = [];

% Peripehral Laminar Recordings in Milo, but in MTC instead 
if (1)
  FolderList{(length(FolderList))+1} = 'Milo_211021';
  FolderList{(length(FolderList))+1} = 'Milo_141021';
  FolderList{(length(FolderList))+1} = 'Milo_121021';
  FolderList{(length(FolderList))+1} = 'Milo_300921';
  FolderList{(length(FolderList))+1} = 'Milo_051021';
  FolderList{(length(FolderList))+1} = 'Milo_071021';
  FolderList{(length(FolderList))+1} = 'Milo_091121';
  FolderList{(length(FolderList))+1} = 'Milo_111121'; 
  FolderList{(length(FolderList))+1} = 'Milo_091221';
  FolderList{(length(FolderList))+1} = 'Milo_141221';
%NEW
   FolderList{(length(FolderList))+1} = 'Milo_180222';
   FolderList{(length(FolderList))+1} = 'Milo_220222';
  FolderList{(length(FolderList))+1} = 'Milo_220322'; 
  FolderList{(length(FolderList))+1} = 'Milo_290322';
  FolderList{(length(FolderList))+1} = 'Milo_010422';
  FolderList{(length(FolderList))+1} = 'Milo_080422';
  FolderList{(length(FolderList))+1} = 'Milo_150422';
  FolderList{(length(FolderList))+1} = 'Milo_190422'; 

end

%%
% Rerun MT files
if (0) % Jude's Laptop
   ProcPath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Processed_MTC';
   StorePath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Info_MTC';
   PlotPath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Plot_Info';
else   % Amy's CPU
    ProcPath = 'C:\Users\abucklaew\Box\MTC\Processed_MT';
    StorePath = 'C:\Users\abucklaew\Box\MTC\MT_WithHartley';
    PlotPath = 'C:\Users\abucklaew\Box\MTC\Plot_Info';
%     ProcPath = 'C:\Users\abucklaew\Box\MTC\Processed_Sprout';
%     StorePath = 'C:\Users\abucklaew\Box\MTC\Info_Sprout';
%     PlotPath = 'C:\Users\abucklaew\Box\MTC\Plot_Info';
end
%******** Folder list for updates
FolderList = [];

% Peripehral Laminar Recordings in Milo, but in MT
if (1)
  FolderList{(length(FolderList))+1} = 'Milo_021221';
  FolderList{(length(FolderList))+1} = 'Milo_040122';
  FolderList{(length(FolderList))+1} = 'Milo_060122';
  FolderList{(length(FolderList))+1} = 'Milo_161221';
  FolderList{(length(FolderList))+1} = 'Milo_301121';
  FolderList{(length(FolderList))+1} = 'Milo_201221';

  FolderList{(length(FolderList))+1} = 'Milo_180122';
  FolderList{(length(FolderList))+1} = 'Milo_280122';
  FolderList{(length(FolderList))+1} = 'Milo_010222';
  FolderList{(length(FolderList))+1} = 'Milo_080222';
  FolderList{(length(FolderList))+1} = 'Milo_250222';
  FolderList{(length(FolderList))+1} = 'Milo_010322';
  FolderList{(length(FolderList))+1} = 'Milo_180322';
  FolderList{(length(FolderList))+1} = 'Milo_050422';
  
%   FolderList{(length(FolderList))+1} = 'Sprout_080922';
%   FolderList{(length(FolderList))+1} = 'Sprout_140922';
%   FolderList{(length(FolderList))+1} = 'Sprout_210922';
%   FolderList{(length(FolderList))+1} = 'Sprout_121022';
%   FolderList{(length(FolderList))+1} = 'Sprout_051022';
%   FolderList{(length(FolderList))+1} = 'Sprout_141022';
%   FolderList{(length(FolderList))+1} = 'Sprout_181022';
%   FolderList{(length(FolderList))+1} = 'Sprout_201022';
%    FolderList{(length(FolderList))+1} = 'Sprout_251022';
end
%% ****** loop over list and do updates
PlotIt = 0;  % if 2, it will save all results as .png
             % if 1, it plots all and accumulates without closing
             % if 0, then no plotting, just save Info struct
for k = 1:length(FolderList)
   FileTag = FolderList{k}; 
   UseFoveal = 0;
   %Analyze_Experiment_SC_Laminar3_MTC(FileTag,ProcPath,StorePath,PlotPath,PlotIt,UseFoveal);
   Analyze_Experiment_AB_Laminar4_MTC(FileTag,ProcPath,StorePath,PlotPath,PlotIt,UseFoveal);
end
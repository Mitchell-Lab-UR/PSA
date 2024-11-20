%%  ********* This is a script to run analysis on multiple experiments
%***
%***
if (1) % Jude's Laptop
   ProcPath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Processed_MT';
   StorePath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\ZInfo_MT';
   PlotPath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Plot_Info';
else   % Amy's CPU
    ProcPath = 'C:\Users\amybu\Box\MTC\Processed_MTC';
    StorePath = 'C:\Users\amybu\Box\MTC\Info_MTC';
    PlotPath = 'C:\Users\amybu\Box\MTC\Plot_Info';
end
%******** Folder list for updates
FolderList = [];


% Peripehral Laminar Recordings in Milo, but in MT
% These are the MT recordings including Hartley stimulus (because there are
% more still beyond just this set)
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
end
%% ****** loop over list and do updates
PlotIt = 0;  % if 2, it will save all results as .png
             % if 1, it plots all and accumulates without closing
             % if 0, then no plotting, just save Info struct
for k = 1:length(FolderList)
   FileTag = FolderList{k}; 
   UseFoveal = 0;
   ZAnalyze_Experiment_Laminar(FileTag,ProcPath,StorePath,PlotPath,PlotIt,UseFoveal);
end
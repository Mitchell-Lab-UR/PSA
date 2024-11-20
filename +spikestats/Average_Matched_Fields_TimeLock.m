function tinfo = Average_Matched_Fields_TimeLock(tinfo, minfo)
% function tinfo = Average_Matched_Fields_TimeLock(tinfo, minfo)
%
%  inputs:  tinfo - larger struct to append a new field
%           minfo - cell of M+1 structs with any set of fields
%                 - it will search all fields, average them into one
%                 - and then append field into tinfo
%           Note: ignore the 1st item of fields (will be empty, or
%                   drawn from an unmatched set)
%  outputs:  tinfo - now with appended struct
%   -- one day I would like this to be general, but for now hack it

   %********
   LockTT = cell(1,2);
   LockUU = cell(1,2);
   LockSU = cell(1,2);
   %********
   AAfit_mu = cell(1,2);
   AAfit_sem = cell(1,2);
   AAtune_mu = cell(1,2);
   AAtune_sem = cell(1,2);
   %********
   Afit_mu = cell(1,2);
   Afit_sem = cell(1,2);
   Atune_mu = cell(1,2);
   Atune_sem = cell(1,2);
   %***************
   Mu = cell(1,2);
   Sem = cell(1,2);
   Num = cell(1,2);
  
   %****************
   M = length(minfo);
   for mm = 1:M
       % mfields = fields(minfo);
       % mfields
       
       for zk = 1:2
           LockTT{zk} = [LockTT{zk} ; minfo{mm}.LockTT{zk}];
           LockUU{zk} = [LockUU{zk} ; minfo{mm}.LockUU{zk}];
           LockSU{zk} = [LockSU{zk} ; minfo{mm}.LockSU{zk}];
           %***** raw fit info
           Mu{zk} = [Mu{zk} ; minfo{mm}.OriTune{zk}.Mu];
           Sem{zk} = [Sem{zk} ; minfo{mm}.OriTune{zk}.Sem];
           Num{zk} = [Num{zk} ; minfo{mm}.OriTune{zk}.Num];
           %********
           AAfit_mu{zk} = [AAfit_mu{zk} ; minfo{mm}.OriTune{zk}.AAfit.mu];
           AAfit_sem{zk} = [AAfit_sem{zk} ; minfo{mm}.OriTune{zk}.AAfit.sem];
           AAtune_mu{zk} = [AAtune_mu{zk} ; minfo{mm}.OriTune{zk}.AAtune.mu];
           AAtune_mu{zk} = [AAtune_mu{zk} ; minfo{mm}.OriTune{zk}.AAtune.sem];
           %*******
           Afit_mu{zk} = [Afit_mu{zk} ; minfo{mm}.OriTune{zk}.Afit.mu];
           Afit_sem{zk} = [Afit_sem{zk} ; minfo{mm}.OriTune{zk}.Afit.sem];
           Atune_mu{zk} = [Atune_mu{zk} ; minfo{mm}.OriTune{zk}.Atune.mu];
           Atune_mu{zk} = [Atune_mu{zk} ; minfo{mm}.OriTune{zk}.Atune.sem];
           %*******       
       end
   end
   
   %******* compute averages
   for zk = 1:2
      %****** 
      tinfo.LockTT{zk} = nanmean(LockTT{zk});
      tinfo.LockUU{zk} = nanmean(LockUU{zk});
      tinfo.LockSU{zk} = nanmean(LockSU{zk});
      %********
      tinfo.OriTune{zk}.Mu = nanmean(Mu{zk});
      tinfo.OriTune{zk}.Sem = nanmean(Sem{zk});
      tinfo.OriTune{zk}.Num = nanmean(Num{zk});
      %********
      tinfo.OriTune{zk}.AAfit.mu = nanmean(AAfit_mu{zk});
      tinfo.OriTune{zk}.AAfit.sem = nanmean(AAfit_sem{zk});
      tinfo.OriTune{zk}.AAtune.mu = nanmean(AAtune_mu{zk});
      tinfo.OriTune{zk}.AAtune.sem = nanmean(AAtune_sem{zk});
      %********
      tinfo.OriTune{zk}.Afit.mu = nanmean(Afit_mu{zk});
      tinfo.OriTune{zk}.Afit.sem = nanmean(Afit_sem{zk});
      tinfo.OriTune{zk}.Atune.mu = nanmean(Atune_mu{zk});
      tinfo.OriTune{zk}.Atune.sem = nanmean(Atune_sem{zk});
      %********
   end
   
end


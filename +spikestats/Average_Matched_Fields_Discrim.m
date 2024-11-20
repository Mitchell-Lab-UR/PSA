function tinfo = Average_Matched_Fields_Discrim(tinfo, minfo)
% function tinfo = Average_Matched_Fields_Discrim(tinfo, minfo)
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
   AucTT = cell(1,2);
   AucUU = cell(1,2);
   AucSU = cell(1,2);
   %********
   DSI_tt = cell(1,2);
   DSI_uu = cell(1,2);
   DSI_su = cell(1,2);
   %********
   Tlocks = length(tinfo.TimeLocks);
   AUC = cell(1,Tlocks);
   AUC_sem = cell(1,Tlocks);
   %****************
   M = length(minfo);
   for mm = 1:M
       % mfields = fields(minfo);
       % mfields
       
       for zk = 1:2
           LockTT{zk} = [LockTT{zk} ; minfo{mm}.LockTT{zk}];
           LockUU{zk} = [LockUU{zk} ; minfo{mm}.LockUU{zk}];
           LockSU{zk} = [LockSU{zk} ; minfo{mm}.LockSU{zk}];
           %******
           AucTT{zk} = [AucTT{zk} ; minfo{mm}.AucTT{zk}];
           AucUU{zk} = [AucUU{zk} ; minfo{mm}.AucUU{zk}];
           AucSU{zk} = [AucSU{zk} ; minfo{mm}.AucSU{zk}];
           %******
           DSI_tt{zk} = [DSI_tt{zk} ; minfo{mm}.DSI_tt{zk}];
           DSI_uu{zk} = [DSI_uu{zk} ; minfo{mm}.DSI_uu{zk}];
           DSI_su{zk} = [DSI_su{zk} ; minfo{mm}.DSI_su{zk}];
       end
       %*******
       for tk = 1:Tlocks
           AUC{tk} = [AUC{tk} ; minfo{mm}.AUC{tk}];
           AUC_sem{tk} = [AUC_sem{tk} ; minfo{mm}.AUC_sem{tk}];
       end
       %*******
       
   end
   
   %******* compute averages
   for zk = 1:2
      %****** 
      tinfo.LockTT{zk} = nanmean(LockTT{zk});
      tinfo.LockUU{zk} = nanmean(LockUU{zk});
      tinfo.LockSU{zk} = nanmean(LockSU{zk});
      %********
      tinfo.AucTT{zk} = nanmean(AucTT{zk});
      tinfo.AucUU{zk} = nanmean(AucUU{zk});
      tinfo.AucSU{zk} = nanmean(AucSU{zk});
      %********
      tinfo.DSI_tt{zk} = nanmean(DSI_tt{zk});
      tinfo.DSI_uu{zk} = nanmean(DSI_uu{zk});
      tinfo.DSI_su{zk} = nanmean(DSI_su{zk});
   end
   %********
   for tk = 1:Tlocks
       tinfo.AUC{1,tk} = nanmean(AUC{tk});
       tinfo.AUC_sem{1,tk} = nanmean(AUC_sem{tk});
   end
   %********
   
end


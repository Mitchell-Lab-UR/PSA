function Info = Compute_SplitDiscrim(Info,Exp,Unit,H,FieldName)
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



% List of trials with Best PFR
       attlist =  Info.AttList;
       apfr = Info.PFR_DotList(attlist);
       mapfr = nanmedian(apfr);
       zz = find( apfr > mapfr);
       Info.AttList_Best = attlist(zz);
       %******
       uttlist = Info.IgnAList;
       upfr = Info.PFR_DotList(uttlist);
       mupfr = nanmedian(upfr);
       zz = find( Info.PFR_DotList(uttlist) > mupfr);
       Info.IgnAList_Best = uttlist(zz);
       %*******
       uttlist = Info.IgnBList;
       upfr = Info.PFR_DotList(uttlist);
       mupfr = nanmedian(upfr);
       zz = find( Info.PFR_DotList(uttlist) > mupfr);
       Info.IgnBList_Best = uttlist(zz);

       
       
% List of trials with worst PFR
       attlist =  Info.AttList;
       apfr = Info.PFR_DotList(attlist);
       mapfr = nanmedian(apfr);
       zz = find( apfr > mapfr);
       Info.AttList_Worst = attlist(zz);
       %******
       uttlist = Info.IgnAList;
       upfr = Info.PFR_DotList(uttlist);
       mupfr = nanmedian(upfr);
       zz = find( Info.PFR_DotList(uttlist) > mupfr);
       Info.IgnAList_Worst = uttlist(zz);
       %*******
       uttlist = Info.IgnBList;
       upfr = Info.PFR_DotList(uttlist);
       mupfr = nanmedian(upfr);
       zz = find( Info.PFR_DotList(uttlist) > mupfr);
       Info.IgnBList_Worst = uttlist(zz);

%%% Then compute AUC for each unit for each condition (Best/Worst PFR)






  
function Info = Compute_SaccadeFix(Info,Exp,H)
%******* 
%******  function Info = Compute_SaccadeFix(Info,Exp,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          H - if integer, produce plots of the results, if figure
%***                  then gives a cell struct with panel handles to plot
%***
%*** Outputs: Info - it updates fields in the Info struct to include 
%***                  information about RT, fix position and drift
%***

%****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0) 
       %****** setup a multi-panel plot (three panels) 
       hh = figure('position',[100 100 1800 600]); 
       HA = subplot('position',[0.150 0.15 0.20 0.7]);
       HB = subplot('position',[0.450 0.15 0.20 0.7]);
       HC = subplot('position',[0.750 0.15 0.20 0.7]);
       %********************
       plot_SacFix_RT(Info,HA);   % plot true vs estimate direction
       plot_SacFix_Pos(Info,HB);  % plot of angle distribution
       plot_SacFix_Vel(Info,HC); % plot of PFR vectors
       if (H == 2)  % store result to png and close
          z = getframe(hh);  % current figure
          uname = [Info.pathplot,filesep,'SacFix_',Info.tagname,'.png'];
          disp(sprintf('Storing image of SacFix graph at %s',uname)); 
          imwrite(z.cdata,uname); % store image for later review
          close(hh);
       end
    else 
       if ~isempty(H)  % H should be a sub-plot handle otherwise
          plot_SacFix_RT(Info,H{1});   % plot true vs estimate direction
          plot_SacFix_Pos(Info,H{2});  % plot of angle distribution
          plot_SacFix_Vel(Info,H{3}); % plot of PFR vectors
       end
    end
    return;
  end

  %*******
  RTlist = cell(1,2);  % reaction times stim onset to saccade 
  FixPos = cell(1,2);  % projection of mean fixation along RF direction, make 2D
  FixVel = cell(1,2);  % projection of drift along RF direction, make 2D
  
  % loop through the list of trials accepted for analyses
  disp('Computing SacFix statistics ...');
  for k = 1:length( Info.Tlist )
    if ~Info.EPhysList(k) 
        continue;   % skip trial if any problem flagging saccade
    end
    disp(sprintf('Processing trial %d of %d',k,length(Info.Tlist)));
    tr = Info.Tlist(k);  % integer index of an included trial
    %**** given trial, grab eye position data for saccade to aperture
    eyeSmo = Exp.D{tr}.eyeSmo;  %smoothed eye data (columns: time, x, y)
    Tonstim = Info.StimOnList(k);    % onset time from Info list
    Tonsac = Info.SacOnList(k);  % offset time from Info list
    if (isnan(Tonstim) || isnan(Tonsac))
        continue;   % skip trial if any problem flagging saccade
    end
    %*** find times in eye trace from stim onset to saccade onset
    ztimes = find( (eyeSmo(:,1) >= Tonstim) & ...
                   (eyeSmo(:,1) <= Tonsac) );
    if isempty(ztimes)
        continue;
    end
    %****** grab the eye traces over the relevant interval
    tt = eyeSmo(ztimes,1);
    xx = eyeSmo(ztimes,2);
    yy = eyeSmo(ztimes,3);
    
    %****** is this towards or away from the RF?
    tg = Info.TargList(k,1); % which of 3 apertures was the target
    targx = Info.TargPosX(k,tg);  % pos in degrees
    targy = Info.TargPosY(k,tg);  % pos in degrees
    rfx = Info.TargPosX(k,1);  % tg == 1 is RF target
    rfy = Info.TargPosY(k,1);  % tg == 1 is RF target
    if (Info.TargList(k,1) == 1)
       astate = 1;
    else
       astate = 2;
    end
    %***** compute the RT *********
    rt = Tonsac - Tonstim;
    RTlist{astate} = [RTlist{astate} ; rt];
    aprj = ( (rfx * mean(xx)) + (rfy * mean(yy)))/norm([rfx,rfy]);
    bprj = ( (rfx * mean(yy)) - (rfy * mean(xx)))/norm([rfx,rfy]);
    FixPos{astate} = [FixPos{astate} ; [aprj,bprj]];
    aprj2 = ( (rfx * median(diff(xx))) + (rfy * median(diff(yy))))/norm([rfx,rfy]);
    aprj2 = aprj2 * 1000;  % degs per sec
    bprj2 = ( (rfx * median(diff(yy))) -(rfy * median(diff(xx))))/norm([rfx,rfy]);
    bprj2 = bprj2 * 1000;  % degs per sec
    FixVel{astate} = [FixVel{astate} ; [aprj2,bprj2]];
    %***********
  end  % loop over the trials
  %****** store results
  Info.RTlist = RTlist;
  %**** negate the median to center two distributions and
  %**** focus on their relative differences ********
  MedPosA = median( [FixPos{1}(:,1) ; FixPos{2}(:,1)]);
  MedPosB = median( [FixPos{1}(:,2) ; FixPos{2}(:,2)]);
  FixPos{1}(:,1) = FixPos{1}(:,1) - MedPosA;
  FixPos{2}(:,1) = FixPos{2}(:,1) - MedPosA;
  FixPos{1}(:,2) = FixPos{1}(:,2) - MedPosB;
  FixPos{2}(:,2) = FixPos{2}(:,2) - MedPosB;
  Info.MedPos = [MedPosA,MedPosB];
  Info.FixPos = FixPos;
  %*****
  MedVelA = median( [FixVel{1}(:,1) ; FixVel{2}(:,1)]);
  MedVelB = median( [FixVel{1}(:,2) ; FixVel{2}(:,2)]);
  FixVel{1}(:,1) = FixVel{1}(:,1) - MedVelA;
  FixVel{2}(:,1) = FixVel{2}(:,1) - MedVelA;
  FixVel{1}(:,2) = FixVel{1}(:,2) - MedVelB;
  FixVel{2}(:,2) = FixVel{2}(:,2) - MedVelB;
  Info.MedVel = [MedVelA,MedVelB];
  Info.FixVel = FixVel;
  %*********
  disp(sprintf('Saccade and fixation statistics finished computing'));
return;

%***** plot reaction times **********
function plot_SacFix_RT(Info,H)
   subplot(H);
   % pfr directions
   attrt = Info.RTlist{1};
   ignrt = Info.RTlist{2};
   vx = 0.08:0.02:0.60;
   atty = hist(attrt,vx);
   itty = hist(ignrt,vx);
   plot(vx,atty,'r.-'); hold on;
   plot(vx,itty,'b.-');
   xlabel('RT (s)');
   ylabel('Counts');
   title(Info.tagname);
   %****
return;

%***** plot mean fixation position **********
function plot_SacFix_Pos(Info,H)
   subplot(H);
   % pfr directions
   att = Info.FixPos{1}(:,1);
   ign = Info.FixPos{2}(:,1);
   vx = -1.0:0.05:1.0;
   atty = hist(att,vx);
   itty = hist(ign,vx);
   plot(vx,atty,'r.-'); hold on;
   plot(vx,itty,'b.-');
   xlabel('Fix Position (degs)');
   ylabel('Counts');
   title(Info.tagname);
   %****
return;

%***** plot mean fixation position **********
function plot_SacFix_Vel(Info,H)
   subplot(H);
   % pfr directions
   att = Info.FixVel{1}(:,1);
   ign = Info.FixVel{2}(:,1);
   vx = -1.5:0.1:1.5;
   atty = hist(att,vx);
   itty = hist(ign,vx);
   plot(vx,atty,'r.-'); hold on;
   plot(vx,itty,'b.-');
   xlabel('Fix Velocity (degs/sec)');
   ylabel('Counts');
   title(Info.tagname);
   %****
   
return;

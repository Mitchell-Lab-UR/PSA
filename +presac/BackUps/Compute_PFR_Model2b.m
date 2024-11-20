function Info = Compute_PFR_Model2b(Info,Exp,H,plotsteps,N)
%******* 
%******  function Info = Compute_PFR_Iqra(Info,Exp,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          H - if integer, produce plots of the results, if figure
%***                  then gives a cell struct with panel handles to plot
%***
%*** Outputs: Info - it updates fields in the Info struct to include 
%***                  information about trial inclusion parameters
%***

%****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0) 
       %****** setup a multi-panel plot (three panels) 
       hh = figure('position',[100 100 1200 400]); 
       HA = subplot('position',[0.10 0.15 0.22 0.70]);
       HB = subplot('position',[0.75 0.15 0.22 0.70]);
       HC = subplot('position',[0.40 0.15 0.25 0.75]);
       plot_PFR_unity(Info,HA,1);   % plot true vs estimate direction, targ
       plot_PFR_unity(Info,HB,2);   % plot true vs estimate direction, other
       plot_PFR_angles(Info,HC);  % plot of angle distribution
      
       %********************
       hh2 = figure('position',[300 300 800 400]); 
       HE = subplot('position',[0.15 0.15 0.70 0.70]);
      
       plot_PFR_traces(Info,HE);  % plot the velocity traces
       if (H == 2)  % store result to png and close
          z = getframe(hh);  % current figure
          uname = [Info.pathplot,filesep,'PFR_',Info.tagname,'.png'];
          disp(sprintf('Storing image of PFR graph at %s',uname)); 
          imwrite(z.cdata,uname); % store image for later review
          close(hh);
          %**********
          z = getframe(hh2);  % current figure
          uname = [Info.pathplot,filesep,'PFR2_',Info.tagname,'.png'];
          disp(sprintf('Storing image of PFR2 graph at %s',uname)); 
          imwrite(z.cdata,uname); % store image for later review
          close(hh2);
       end
    else 
       if ~isempty(H)  % H should be a sub-plot handle otherwise
          plot_PFR_unity(Info,H{1},1);   % plot true vs estimate direction
          plot_PFR_unity(Info,H{2},2);   % plot true vs estimate direction
          plot_PFR_angles(Info,H{3});  % plot of angle distribution
          plot_PFR_traces(Info,H{4});  % plot the velocity traces
       end
    end
    return;
  end

  % define parameters for doing the PFR analysis
  VELTHRESH = 40;  % threshold to flag second saccades
  CURVETHRESH = 1.2;  % threshold to flag curved saccades
  %********
  OpenLoopStart = 0.02;
  OpenLoopEnd = 0.10;
  PostStart = 0.00;
  PostEnd = 0.20;
  PreStart = -0.10;
  PreEnd = 0.00;
  StimSpeed = 15;  % would be better to read from Exp file
  Max_Eye_Speed = 20; % if great, NaN from smooth traces
  Max_Eye_Gap = 10;  % NaN +/- this many samples from over Max_Eye
  GSIG = 5;  % smooth the velocity profiles extra
  %*****
  Info.pfr_rt = [];
  NT = Info.NumTargs;
  Info.vec = cell(1,NT); % vectors by attended location
  Info.svec =  cell(1,NT); % vectors by spatial location
  Info.avec = cell(1,NT);
  Info.openvec = []; 
  Info.pfrs = cell(1,NT);  % for each target location, 3 of them
                             % it will be a T x 4 matrix, where
                             %  1st column is dot product on target motion
                             %  2nd column is cross product
                             %  3rd column gives the target used (by loc)
                             %  4th column gives the direction of motion
                             %                    (at the loc)
                             % 1st matrix is for the saccade target
                             % 2nd matrix is the first distractor
                             %     (clockwise of the saccade target)
                             % 3rd matrix is the second distractor
  %**** Store the parameters used in analysis                           
  Info.OpenLoopStart = OpenLoopStart;
  Info.OpenLoopEnd = OpenLoopEnd;
  Info.OpenLoopSpeed = StimSpeed;
  Info.PostStart = PostStart;
  Info.PostEnd = PostEnd;
  Info.PreStart = PreStart;
  Info.PreEnd = PreEnd;
  %***********
  % plotsteps = 0;  % plot trial by trial to confirm it is working
  %***********
  
  %***** NORMED:  assume each target location causes a default drift
  %************:  run through all trials, store PFR by location and motion
  %************:  compute the weighted average drift at each location
  Info.LocationBias = cell(1,NT);  % each cell a list of x,y vec PFR
  Info.MeanLocBias = cell(1,NT);
  for k = 1:length( Info.Tlist )
    if ~Info.EPhysList(k) 
        continue;   % skip trial if any problem flagging saccade
    end
    disp(sprintf('Processing Location Bias: trial %d of %d',k,length(Info.Tlist)));
    tr = Info.Tlist(k);  % integer index of an included trial
    %**** given trial, grab eye position data for saccade to aperture    
    eyeSmo = Exp.D{tr}.eyeSmo;  %smoothed eye data (columns: time, x, y)
    Tonsac = Info.SacOnList(k);    % onset time from Info list
    Toffsac = Info.SacOffList(k);  % offset time from Info list
    Tonstim = Info.StimOnList(k);    % onset time from Info list
    if (isnan(Tonsac) || isnan(Toffsac))
        continue;   % skip trial if any problem flagging saccade
    end
    
    %compute PFR vector
    openloop = find((eyeSmo(:,1) >= (Toffsac+OpenLoopStart)) & (eyeSmo(:,1) <= (Toffsac+OpenLoopEnd)));
    ott = eyeSmo(openloop,1);
    oxx = eyeSmo(openloop,2);
    oyy = eyeSmo(openloop,3);
    
    %compute startposition and open loop vector (openvec)
    if isempty(openloop)
        continue;
    end
    startxx = eyeSmo(openloop(1),2);
    startyy = eyeSmo(openloop(1),3);
    startposition = [startxx, startyy];
    pxx = (eyeSmo(openloop(end),2) - startxx); % x position difference
    pyy = (eyeSmo(openloop(end),3) - startyy); % y position difference
    
    %******* store vectors, is there a net drift per location
    tg = Info.TargList(k,1); % which of 3 apertures was the target 
    Info.LocationBias{tg} = [Info.LocationBias{tg} ; [pxx,pyy] ];
  end
  for k = 1:NT
      Info.MeanLocBias{k} = median(Info.LocationBias{k});  % use median, reduce outlier
  end
  %**************** 
 
  if (plotsteps == 1)
            %**********plot some of the data
            h = figure(1); 
            set(h,'position',[100 100 1200 800]);  hold off;
  end
  
  %***** store velocity traces
  Info.VelTraces = cell(2,2);
  % loop through the list of trials accepted for analyses
  disp('Computing PFR statistics ...');
  Info.PFR_DotList = nan(size(Info.Tlist));  % store dot prods 
  Info.TrialFlag = [];
  for k = 1:length( Info.Tlist )
    if ~Info.EPhysList(k) 
        continue;   % skip trial if any problem flagging saccade
    end
    disp(sprintf('Processing trial %d of %d',k,length(Info.Tlist)));
    tr = Info.Tlist(k);  % integer index of an included trial
    %**** given trial, grab eye position data for saccade to aperture    
    eyeSmo = Exp.D{tr}.eyeSmo;  %smoothed eye data (columns: time, x, y)
    Tonsac = Info.SacOnList(k);    % onset time from Info list
    Toffsac = Info.SacOffList(k);  % offset time from Info list
    Tonstim = Info.StimOnList(k);    % onset time from Info list
    if (isnan(Tonsac) || isnan(Toffsac))
        continue;   % skip trial if any problem flagging saccade
    end
    %*** find times in eye trace -0.1 to +0.1 around saccade offset
    ztimes = find( (eyeSmo(:,1) >= (Toffsac-0.10)) & ...
                   (eyeSmo(:,1) <= (Toffsac+0.12)) );
    if isempty(ztimes)
        continue;
    end
    %****** grab the eye traces over the relevant interval
    tt = eyeSmo(ztimes,1);
    xx = eyeSmo(ztimes,2);
    yy = eyeSmo(ztimes,3);
    
    %compute PFR vector
    openloop = find((eyeSmo(:,1) >= (Toffsac+OpenLoopStart)) & (eyeSmo(:,1) <= (Toffsac+OpenLoopEnd)));
    ott = eyeSmo(openloop,1);
    oxx = eyeSmo(openloop,2);
    oyy = eyeSmo(openloop,3);
    
    %compute startposition and open loop vector (openvec)
    if isempty(openloop)
        continue;
    end
    startxx = eyeSmo(openloop(1),2);
    startyy = eyeSmo(openloop(1),3);
    startposition = [startxx, startyy];    
    pxx = (eyeSmo(openloop(end),2) - startxx); % x position difference
    pyy = (eyeSmo(openloop(end),3) - startyy); % y position difference
    % openvec = [pxx,pyy];
   
    %******* determine mean drift per target location
    tg = Info.TargList(k,1); % which of 3 apertures was the target
    mxx = Info.MeanLocBias{tg}(1);
    myy = Info.MeanLocBias{tg}(2);  % median PFR shift per location
    vmxx = mxx/(OpenLoopEnd-OpenLoopStart);
    vmyy = myy/(OpenLoopEnd-OpenLoopStart); % median vel vector per loc
    %*********
    openvec = [(pxx-mxx),(pyy-myy)];
    
    %********** compute absolute velocity over the same window (ztimes) as
    %********** the position traces were computed
    abvxx = [];
    abvyy = [];
    abvsp = [];
    abvtt = [];
    for t = 2:length(xx)
       velx = (xx(t) - xx(t-1))/(tt(t)-tt(t-1));  % dva / secs
       vely = (yy(t) - yy(t-1))/(tt(t)-tt(t-1));  % dva / secs
       vabs = norm([velx,vely]);
       abvtt = [abvtt ; ((tt(t)+tt(t-1))/2)];  % mid-point of time between the two points used to compute velocity
       abvxx = [abvxx velx];
       abvyy = [abvyy vely];
       abvsp = [abvsp vabs];
    end
    %******* Apply additional smoothing for velocity
    abvxx = spikestats.gauss_smooth(abvxx,GSIG);
    abvyy = spikestats.gauss_smooth(abvyy,GSIG);
    abvsp = spikestats.gauss_smooth(abvsp,GSIG);
   
    %********** compute interval to analyze velocity traces
    ploop = find((eyeSmo(:,1) >= (Toffsac+PostStart)) & (eyeSmo(:,1) <= (Toffsac+PostEnd)));
    plop = ploop(1) + (0:floor((PostEnd-PostStart)/0.001));
    if (max(plop) > size(eyeSmo,1) )
       plxx = NaN(1,length(plop));
       plyy = NaN(1,length(plop));
       pltt = 1:length(plop);
       disp('Failed to reach limit of eye position');
    else    
       plxx = eyeSmo(plop,2);
       plyy = eyeSmo(plop,3);
       pltt = eyeSmo(plop,1);
    end
    %***** compute the velocity traces
    vxx = [];
    vyy = [];
    vsp = [];
    vtt = [];
    for t = 2:length(plxx)
       velx = (plxx(t) - plxx(t-1))/(pltt(t)-pltt(t-1));  % dva / secs
       vely = (plyy(t) - plyy(t-1))/(pltt(t)-pltt(t-1));  % dva / secs
       vabs = norm([velx,vely]);
       vtt = [vtt ; ((pltt(t)+pltt(t-1))/2)];  % mid-point of time between the two points used to compute velocity
       vxx = [vxx velx];
       vyy = [vyy vely];
       vsp = [vsp vabs];
    end
    %******* Apply additional smoothing for velocity
    vxx = spikestats.gauss_smooth(vxx,GSIG);
    vyy = spikestats.gauss_smooth(vyy,GSIG);
    vsp = spikestats.gauss_smooth(vsp,GSIG);
    
    %********** compute interval to analyze velocity traces
    ploop2 = find((eyeSmo(:,1) >= (Tonsac+PreStart)) & (eyeSmo(:,1) <= (Tonsac+PreEnd)));
    plop2 = ploop2(1) + (0:floor((PreEnd-PreStart)/0.001));
    if (max(plop2) > size(eyeSmo,1) )
       plxx2 = NaN(1,length(plop2));
       plyy2 = NaN(1,length(plop2));
       pltt2 = 1:length(plop2);
       disp('Failed to reach limit of eye position');
    else    
       plxx2 = eyeSmo(plop2,2);
       plyy2 = eyeSmo(plop2,3);
       pltt2 = eyeSmo(plop2,1);
    end
    %***** compute the velocity traces
    vxx2 = [];
    vyy2 = [];
    vsp2 = [];
    for t = 2:length(plxx2)
       velx2 = (plxx2(t) - plxx2(t-1))/(pltt2(t)-pltt2(t-1));  % dva / secs
       vely2 = (plyy2(t) - plyy2(t-1))/(pltt2(t)-pltt2(t-1));  % dva / secs
       vabs2 = norm([velx2,vely2]);
       vxx2 = [vxx2 velx2];
       vyy2 = [vyy2 vely2];
       vsp2 = [vsp2 vabs2];
    end
    %******* Apply additional smoothing for velocity
    vxx2 = spikestats.gauss_smooth(vxx2,GSIG);
    vyy2 = spikestats.gauss_smooth(vyy2,GSIG);
    vsp2 = spikestats.gauss_smooth(vsp2,GSIG);
   
    %***** go through and NaN away velocities above max speed
    zvxx = [vxx2 NaN vxx];
    zvyy = [vyy2 NaN vyy];
    zvsp = [vsp2 NaN vsp];  % absoluted velocity
    %**********
    nvxx = zvxx;
    nvyy = zvyy;
    nvsp = zvsp;
    for t = 1:length(nvsp)
        if ~isnan(nvsp(t)) && (nvsp(t) > Max_Eye_Speed)
           ta = max(1,(t-Max_Eye_Gap));
           tb = length(nvsp); % min(length(nvsp),(t+Max_Eye_Gap));
           nvxx(ta:tb) = NaN; % throw out rest of trial if saccade
           nvyy(ta:tb) = NaN;
        end
    end
    %*******
    
    %******** Get the vectors by spatial location (not by target)  
    motion1 = Info.TargMotion(k,1);  % motion in degrees for first spatial target
    motion2 = Info.TargMotion(k,2);  % motion in degrees (second)
    motion3 = Info.TargMotion(k,3);  % motion in degrees (third)

    %*********
    endspeed = StimSpeed * (OpenLoopEnd-OpenLoopStart);
    sango = motion1 * pi/180;
    svec = endspeed * [cos(sango), sin(sango)];
    sango2 = motion2 * pi/180;
    svec2 = endspeed * [cos(sango2), sin(sango2)];
    sango3 = motion3 * pi/180;
    svec3 = endspeed * [cos(sango3), sin(sango3)];
    %***********
    %******* store information into svec matrices

    %******* basic RF information per location, here provide non-targ vecs
    %*******  but zero the target vector (will be non-zero in other param)
    tg = Info.TargList(k,1);
    if (tg == 1)
      Info.vec{1} = [Info.vec{1} ; [0,0]];
    else
      Info.vec{1} = [Info.vec{1} ; svec ];
    end
    if (tg == 2)
      Info.vec{2} = [Info.vec{2} ; [0,0]];      
    else
      Info.vec{2} = [Info.vec{2} ; svec2];
    end
    if (tg == 3)
      if (NT > 2)     
          Info.vec{3} = [Info.vec{3} ; [0,0]];
      end
    else
      if (NT > 2)     
          Info.vec{3} = [Info.vec{3} ; svec3];
      end      
    end
    %****** for saccade target, give vector, else zero it
    if (tg == 1)
       Info.svec{1} = [Info.svec{1} ; svec];
       Info.svec{2} = [Info.svec{2} ; [0,0]];
       if (NT > 2)
           Info.svec{3} = [Info.svec{3} ; [0,0]];
       end
    end
    if (tg == 2)
       Info.svec{1} = [Info.svec{1} ; [0,0]];
       Info.svec{2} = [Info.svec{2} ; svec2];
       if (NT > 2)
           Info.svec{3} = [Info.svec{3} ; [0,0]];
       end    
    end
    if (tg == 3)
       Info.svec{1} = [Info.svec{1} ; [0,0]];
       Info.svec{2} = [Info.svec{2} ; [0,0]];
       if (NT > 2)
           Info.svec{3} = [Info.svec{3} ; svec3];
       end          
    end
    %***************

    %******** Get the target motions
    tg = Info.TargList(k,1); % which of 3 apertures was the target
    tgmotion = Info.TargMotion(k,tg);  % motion in degrees
    if (NT == 2)
      if (tg == 2)
          tg2 = 1;
      else
          tg2 = 2;
      end
      tgmotion2 = Info.TargMotion(k,tg2);  % motion in degrees
      tgmotion3 = 0;  % fix to non-informative
    else  
        tg2 = tg+1;
        if (tg2 > 3)
            tg2 = tg2 - 3;
        end
        tgmotion2 = Info.TargMotion(k,tg2);  % motion in degrees
        tg3 = tg+2;
        if (tg3 > 3)
            tg3 = tg3 - 3;
        end
        tgmotion3 = Info.TargMotion(k,tg3);
    end
    %*********
    endspeed = StimSpeed * (OpenLoopEnd-OpenLoopStart);
    ango = tgmotion * pi/180;
    tgvec = endspeed * [cos(ango), sin(ango)];
    ango2 = tgmotion2 * pi/180;
    tgvec2 = endspeed * [cos(ango2), sin(ango2)];
    ango3 = tgmotion3 * pi/180;
    tgvec3 = endspeed * [cos(ango3), sin(ango3)];
    %***********
    
    % compute PFR in the open loop secs after offset
    openlength = sqrt( openvec(1)^2 + openvec(2)^2 );  % size of PFR vec
    tglength = sqrt( tgvec(1)^2 + tgvec(2)^2);
    dot = ((tgvec(1)*openvec(1)) + (tgvec(2)*openvec(2)))/tglength;
    cro = -((tgvec(1)*openvec(2)) - (tgvec(2)*openvec(1)))/tglength;
    tglength2 = sqrt( tgvec2(1)^2 + tgvec2(2)^2);
    dot2 = ((tgvec2(1)*openvec(1)) + (tgvec2(2)*openvec(2)))/tglength;
    cro2 = -((tgvec2(1)*openvec(2)) - (tgvec2(2)*openvec(1)))/tglength;
    tglength3 = sqrt( tgvec3(1)^2 + tgvec3(2)^2);
    dot3 = ((tgvec3(1)*openvec(1)) + (tgvec3(2)*openvec(2)))/tglength;
    cro3 = -((tgvec3(1)*openvec(2)) - (tgvec3(2)*openvec(1)))/tglength;
    
    %*******
    Info.PFR_DotList(k) = dot;  % store dot product of PFR per trial
    
    %******** project velocity onto motion vector, to get vdot
    tglength = sqrt( tgvec(1)^2 + tgvec(2)^2);  % target motion vector
    vdot = ((tgvec(1) * vxx) + (tgvec(2) * vyy))/(tglength * StimSpeed);
    nvdot = ((tgvec(1) * nvxx) + (tgvec(2) * nvyy))/(tglength * StimSpeed);
    %********
    if isfield(Exp.D{tr}.PR,'targori')
         if isnan( Exp.D{tr}.PR.changori )  % if NaN, it blanked
               Info.VelTraces{1,2} = [Info.VelTraces{1,2} ; vdot];
               Info.VelTraces{2,2} = [Info.VelTraces{2,2} ; nvdot];
         else
            if (Exp.D{tr}.PR.changori == Exp.D{tr}.PR.targori) %here stayed the same
               Info.VelTraces{1,1} = [Info.VelTraces{1,1} ; vdot];
               Info.VelTraces{2,1} = [Info.VelTraces{2,1} ; nvdot];
            end
         end
    end
    %*******
    
    rt = Tonsac - Tonstim;
    Info.pfr_rt = [Info.pfr_rt ; rt];
    %******* store information into prfs matrices
    % Info.pfrs = cell{1,3};  % for each target location, 3 of them
    Info.pfrs{1} = [Info.pfrs{1} ; [dot,cro,tg,tgmotion]];
    Info.pfrs{2} = [Info.pfrs{2} ; [dot2,cro2,tg2,tgmotion2]];
    if (NT > 2)
       Info.pfrs{3} = [Info.pfrs{3} ; [dot3,cro3,tg3,tgmotion3]];
    end
    
    %******* store information into vec matrices
    
    %**** See how I recoded this above on lines 325...
    Info.avec{1} = [Info.avec{1} ; tgvec ];
    Info.avec{2} = [Info.avec{2} ; tgvec2];
    if (NT > 2)
       Info.avec{3} = [Info.avec{3} ; tgvec3];
    end
    Info.openvec = [Info.openvec ; openvec];
    
    %******* Find a second based on velocity in the post-saccade period
    secondsac = 0;  % this would equal 1 if it happened
    zpostsac = find( (abvtt >= (Toffsac+OpenLoopStart)) & (abvtt < (Toffsac+OpenLoopEnd)));
    maxvel = max(abvsp(zpostsac));  % max vel in the open loop period
    if (maxvel > VELTHRESH)
        secondsac = 1;
    end
    
    %****** Find if a main saccade is curved too much from straight
    failedsac = 0;
    curvesac = 0;
    %**** first get straight line distance
    %*** find the start position
    zstart = find( (tt >= Tonsac) & (tt < (Tonsac+0.01)) );
    %*** find the end position
    zend = find( (tt >= Toffsac) & (tt < (Toffsac+0.01)) );
    if isempty(zstart) | isempty(zend)
        failedsac = 1;
    else
      xstart = xx(zstart(1));
      ystart = yy(zstart(1));
      xend = xx(zend(1));
      yend = yy(zend(1));
      %**** compute distance of a straight line
      sacdist = norm([(xend-xstart),(yend-ystart)]);
      %**** now integrated the path of saccade for its length
      pathsum = 0;
      for tk = zstart(2):zend(1)
        dist = norm([(xx(tk)-xx(tk-1)),(yy(tk)-yy(tk-1))]);
        pathsum = pathsum + dist;
      end
      curveportion = (pathsum/sacdist);
      if (curveportion > CURVETHRESH)
        curvesac = 1;
      end
    end
    
    %******* Store our flagging of trials
    trialflag = 1;  % default is a good trial
    if (secondsac == 1)
        trialflag = 2;
    end
    if (curvesac == 1)
        trialflag = 3;
    end
    if (failedsac == 1)
        trialflag = 4;
    end
    Info.TrialFlag(k) = trialflag;
    %************************
    
    %****** LAST is to plot the results on a trial by trial to debug
    if (plotsteps == 1)
            %**********plot some of the data
            % h = figure(1); 
            % set(h,'position',[100 100 800 800]);  hold off;
            
            %**** plot the spatial eye trace from fixation
            subplot(2,3,1); hold off;
            plot(xx,yy,'k.-'); hold on;
            h = plot(oxx,oyy,'k.-');
            if (secondsac == 1)   % second saccade happened
                 set(h,'Color',[1.0,1.0,0.0]);  % label open loop trace in gray
            else
                if (curvesac == 1)
                   set(h,'Color',[1.0,0.5,0.0]);  % label open loop trace in gray       
                else
                   set(h,'Color',[1.0,0.0,0.0]);  % label open loop trace in gray
                end
            end  
            set(h,'Markersize',10);  % label open loop trace in gray
            plot(0,0,'r+');  % fixation point
            wid = 8;
            axis([-wid wid -wid wid]);  % bound in -10 to +10 window
            xlabel('X Eye Position');
            ylabel('Y Eye Position');
            title(sprintf('Spatial Plots:Trial(%d,%d)',tr,k));
            %****** Inserting new information here about targets positions and
            %****** and their motion directions
            for ti = 1:size(Info.TargMotion,2)
              tx = Info.TargPosX(k,ti);
              ty = Info.TargPosY(k,ti);
              modir = ( Info.TargMotion(k,ti) * (pi/180));
              plot(tx,ty,'bo'); 
              h = plot([tx,tx+cos(modir)],[ty,ty+sin(modir)],'b:');
              set(h,'Linewidth',2);
            end
            wid = 2;
            xa = Info.TargPosX(k,tg);
            ya = Info.TargPosY(k,tg);
            plot([(xa-wid),(xa-wid),(xa+wid),(xa+wid),(xa-wid)],...
                 [(ya-wid),(ya+wid),(ya+wid),(ya-wid),(ya-wid)],'g-');
            %************** 
            
            %**** plot the time-course of traces 
            subplot(2,3,3); hold off;
            plot(tt,xx,'r.-'); hold on;
            plot(tt,yy,'b.-');
            axis([tt(1) tt(end) -10 10]);
            plot([Tonsac,Tonsac],[-10,10],'k-');
            plot([Toffsac,Toffsac],[-10,10],'k-');
            xlabel('Time (secs)');
            ylabel('Position');
            title('Traces over time (H-red,V-blue)');
            
            %**** plot the time-course of traces 
            subplot(2,3,6); hold off;
            plot(abvtt,abvsp,'k-'); hold on;
            axis([abvtt(1) abvtt(end) -10 500]);
            plot([Tonsac,Tonsac],[-10,500],'k-');
            plot([Toffsac,Toffsac],[-10,500],'k-');
            xlabel('Time (secs)');
            ylabel('Velocity (dva/s)');
            
            %******* show a zoom in of the traces
            subplot(2,3,2); hold off;
            %******* replots what was in the first panel
            if (1)
              plot(xx,yy,'k.-'); hold on; 
              h = plot(oxx,oyy,'k.-');
              if (secondsac == 1)   % second saccade happened
                 set(h,'Color',[1.0,1.0,0.0]);  % label open loop trace in gray
              else
                if (curvesac == 1)
                   set(h,'Color',[1.0,0.5,0.0]);  % label open loop trace in gray       
                else
                   set(h,'Color',[1.0,0.0,0.0]);  % label open loop trace in gray
                end
              end  
              set(h,'Markersize',10);
              plot(0,0,'r+');  % fixation point
              xlabel('X Eye Position');
              ylabel('Y Eye Position');
              title('Zoom in at Target');
              %****** Inserting new information here about targets positions and
              %****** and their motion directions
              for ti = 1:size(Info.TargMotion,2)
                tx = Info.TargPosX(k,ti);
                ty = Info.TargPosY(k,ti);
                modir = ( Info.TargMotion(k,ti) * (pi/180));
                plot(tx,ty,'bo'); 
                h = plot([tx,tx+cos(modir)],[ty,ty+sin(modir)],'b:');
                set(h,'Linewidth',2);
              end
            end
            %***** but zooms in on the saccade target end-point
            wid = 2;
            xa = Info.TargPosX(k,tg);
            ya = Info.TargPosY(k,tg);
            axis([xa-wid xa+wid ya-wid ya+wid]);  % bound in -10 to +10 window
            wid = 1.95;
            plot([(xa-wid),(xa-wid),(xa+wid),(xa+wid),(xa-wid)],...
                 [(ya-wid),(ya+wid),(ya+wid),(ya-wid),(ya-wid)],'g-');
            title(sprintf('2ndsac(%d) Curve(%d)',secondsac,curvesac));
            
            %plot the target motion and the PFR vector with same origin to compare
            subplot(2,3,4); hold off;
            plot(0,0,'ko'); hold on;
            plot([0,tgvec(1)],[0,tgvec(2)],'b:'); 
            plot([0,openvec(1)],[0,openvec(2)],'r-');
            axis(1.2 * [-endspeed endspeed -endspeed endspeed]);
            xlabel('x velocity');
            ylabel('y velocity');
            title('PFR and target motion vectors');

            % aligned the PFR vector along with the motion being vertical 
            subplot(2,3,5); % hold off;
            plot(0,0,'ko'); hold on;
            plot([0,0],[0,tglength],'b-'); %,'Linewidth',2);
            plot([0,cro],[0,dot],'r-');
            axis(1.2*[-endspeed endspeed -endspeed endspeed]);
            xlabel('x velocity');
            ylabel('y velocity');
            title('PFR vector aligned to motion');
        
            
        %***** stop to check each trial
        input('Press enter to continue');
    end  % if plotsteps    
  end  % loop over the trials

  
  %****** as a final session stat, compute median angular error for target
  %****** and for the non-targets per session
  dots = Info.pfrs{1}(:,1); 
  cros = Info.pfrs{1}(:,2); 
  tglength =  Info.OpenLoopSpeed * (Info.OpenLoopEnd-Info.OpenLoopStart);
  mago = abs( dots + i * cros);
  zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
  angos = angle( dots(zz) +  i * cros(zz) ) * (180/pi);  % go from radians to degs
  this_angos = angos;
  
  Info.PFR_TargErr = median(abs(angos));
  %*** resultant vector
  udot = mean(dots(zz));
  ucro = mean(cros(zz));
  Info.PFR_TargResV = [udot,ucro];
  %**** compute histograms of PFR error
  vx = -180:10:180;
  vy = hist(angos,vx);
  vy = vy / sum(vy);   % normalize
  %*******
  Info.PFR_TargHistX = vx;
  Info.PFR_TargHist = vy;

  %*******
  if (NT > 2)
    dots = [ Info.pfrs{2}(:,1); Info.pfrs{3}(:,1)];
    cros = [ Info.pfrs{2}(:,2); Info.pfrs{3}(:,2)];
  else
    dots = [ Info.pfrs{2}(:,1)];
    cros = [ Info.pfrs{2}(:,2)];    
  end
  mago = abs( dots + i * cros);
  zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
  angos = angle( dots(zz) +  i * cros(zz) ) * (180/pi);  % go from radians to degs
  Info.PFR_DistErr = median(abs(angos));
  vy = hist(angos,vx);
  vy = vy / sum(vy);   % normalize
  %*** resultant vector
  udot = mean(dots(zz));
  ucro = mean(cros(zz));
  Info.PFR_DistResV = [udot,ucro];
  %*******
  Info.PFR_DistHist = vy;
  
  %******* Velocity traces to store for pooling *************
  Info.VelTraces_U1= nanmean(Info.VelTraces{2,1});
  Info.VelTraces_S1 = nanstd(Info.VelTraces{2,1})/sqrt(size(Info.VelTraces{2,1},1));
  Info.VelTraces_U2 = nanmean(Info.VelTraces{2,2});
  Info.VelTraces_S2 = nanstd(Info.VelTraces{2,2})/sqrt(size(Info.VelTraces{2,2},1));
 
%******** last step is to plot results if desired
disp(sprintf('PFR statistics finished computing: Targ(%5.2f) Ign(%5.2f)',...
      Info.PFR_TargErr,Info.PFR_DistErr)); 
         
 %******** Try to use lm to fit PFR trial by trial here ********
  
  pfr_x =  Info.openvec(:,1); % PFR vector
  pfr_y =  Info.openvec(:,2);
   
  targ_x = Info.vec{1}(:,1); % Attended Target vectors 
  targ_y = Info.vec{1}(:,2); 
  other1_x = Info.vec{2}(:,1); 
  other1_y = Info.vec{2}(:,2); 
  other2_x = Info.vec{3}(:,1); 
  other2_y = Info.vec{3}(:,2); 
  
  loc1_x = Info.svec{1}(:,1); % Spatial Target vectors 
  loc1_y = Info.svec{1}(:,2); 
  loc2_x = Info.svec{2}(:,1); 
  loc2_y = Info.svec{2}(:,2); 
  loc3_x = Info.svec{3}(:,1); 
  loc3_y = Info.svec{3}(:,2); 
  
  %****** try to debug on old model
  ztarg_x = Info.avec{1}(:,1); % Attended Target vectors 
  ztarg_y = Info.avec{1}(:,2); 
  zother1_x = Info.avec{2}(:,1); 
  zother1_y = Info.avec{2}(:,2); 
  zother2_x = Info.avec{3}(:,1); 
  zother2_y = Info.avec{3}(:,2); 
  
  
  %features_x = [targ_x,other1_x,other2_x,loc1_x,loc2_x,loc3_x];
  %features_y = [targ_y,other1_y,other2_y,loc1_y,loc2_y,loc3_y];
 
  features_x = [ztarg_x,zother1_x,zother2_x];
  features_y = [ztarg_y,zother1_y,zother2_y];
 
  mdl_x = fitlm(features_x,pfr_x); % lm model for x vector
  mdl_y = fitlm(features_y,pfr_y); % lm model for y vector
  
  p_x = predict(mdl_x,features_x); % predict new model directions for x vectors
  p_y = predict(mdl_y,features_y); % predict new model directions for y vectors
  
  mdl_x
  mdl_y
  
  %******* force same weight values for x,y (might be n
  if (0)
      beta_x = mdl_x.Coefficients.Estimate;
      beta_y = mdl_y.Coefficients.Estimate;
      ux = mean(beta_x(2:4));
      uy = mean(beta_y(2:4));
      beta_u = 0.5*(beta_x + beta_y);
      beta_x(2:4) = beta_u(2:4);
      beta_y(2:4) = beta_u(2:4);
      % beta_x(2:7) = beta_u(2:7) * (ux / 0.5*(ux + uy));
      % beta_y(2:7) = beta_u(2:7) * (uy / 0.5*(ux + uy));
      %**********
      myN = size(features_x,1);
      augfeatures_x = [ones(myN,1) features_x];
      augfeatures_y = [ones(myN,1) features_y];
      %*******
      p_x = augfeatures_x * beta_x;
      p_y = augfeatures_y * beta_y;
      %*******
      mdl_x
      mdl_y
      beta_u
  end
  
%   noise_x = pfr_x-p_x; % find the residuals for each vector from true pfr
%   noise_y = pfr_y-p_y;
%    
%   for k = 1:N; 
%       New_p_x(:,k) = p_x + noise_x(randperm(length(noise_x))); % randomize residual and add them back to predicted vectors
%       New_p_y(:,k) = p_y + noise_y(randperm(length(noise_y)));
%   end 
%       
%   for k = 1:N; 
%       p_angos(:,k) = wrapTo180(atan2d(New_p_y(:,k),New_p_x(:,k)) - atan2d(targ_y,targ_x)); % recalcuate angular error, but need to subtract out from targ vector
%   end
%   

  p_angos = wrapTo180(atan2d(p_y,p_x) - atan2d(targ_y,targ_x)); % recalcuate angular error, but need to subtract out from targ vector
  
  Info.PFR_Model_TargErr = median(median(abs(p_angos))); % Save median error for the model
  
  %**** compute histograms of PFR error
  vx = -180:10:180;
  vy = hist(p_angos,vx);
  vy = vy / sum(vy);   % normalize
  %*******
  Info.PFR_Model_TargHist = vy;  
%   
%   for k =1:N;
%       p_angos1(:,k) = wrapTo180(atan2d(New_p_y(:,k),New_p_x(:,k))- atan2d(other1_y,other1_x));
%       p_angos2(:,k) = wrapTo180(atan2d(New_p_y(:,k),New_p_x(:,k))- atan2d(other2_y,other2_x));
%   end
  
      p_angos1 = wrapTo180(atan2d(p_y,p_x)- atan2d(other1_y,other1_x));
      p_angos2 = wrapTo180(atan2d(p_y,p_x)- atan2d(other2_y,other2_x));
  
  
  Info.PFR_Model_DistErr = median(median(abs([p_angos1; p_angos2])));
  
  Info.p_angos = p_angos(:,1);  %% only save out one model run? 
  Info.p_angos1 = p_angos1(:,1);
  Info.p_angos2 = p_angos2(:,1); 
  
  %**** compute histograms of PFR error
   vx = -180:10:180;
   vy = hist([p_angos1; p_angos2],vx);
   vy = vy / sum(vy);   % normalize
  %*******
  Info.PFR_Model_DistHist = vy;
     
return;








%***** plot the true motion versus estimated direction on a unity line          
function plot_PFR_unity(Info,H,Targ)
  
   subplot(H);
   % pfr directions
   zz = find( Info.TrialFlag == 1);
   %*********
   if (Targ == 1)
     dots = Info.pfrs{1}(:,1); 
     cros = Info.pfrs{1}(:,2); 
     tangs = Info.pfrs{1}(:,4); % actual motion of target apertures
   else
     NT = Info.NumTargs;  
     if (NT > 2)
       dots = [ Info.pfrs{2}(:,1); Info.pfrs{3}(:,1)];
       cros = [ Info.pfrs{2}(:,2); Info.pfrs{3}(:,2)];
       tangs = [Info.pfrs{2}(:,4); Info.pfrs{3}(:,4)]; % actual motion of target apertures
     else
       dots = [ Info.pfrs{2}(:,1)];
       cros = [ Info.pfrs{2}(:,2)];    
       tangs = [Info.pfrs{2}(:,4)]; % actual motion of target apertures
     end
   end
   %********
   angos = angle( dots +  i * cros ) *  (180/pi);  % go from radians to degs
   tangs = tangs; % match vectors
   
   p_angos = Info.p_angos + Info.pfrs{1}(:,4);
   p_angos1 = Info.p_angos1 + Info.pfrs{2}(:,4);
   p_angos2 = Info.p_angos2 + Info.pfrs{3}(:,4);
   
      zz2 = find(p_angos > 360);  % wrap around circle)
   if ~isempty(zz2)
       p_angos(zz2) = p_angos(zz2) - 360;
   end
   zz2 = find(p_angos < 0);  % wrap around circle)
   if ~isempty(zz2)
       p_angos(zz2) = p_angos(zz2) + 360;
   end
   
      zz2 = find(p_angos1 > 360);  % wrap around circle)
   if ~isempty(zz2)
       p_angos1(zz2) = p_angos1(zz2) - 360;
   end
   zz2 = find(p_angos1 < 0);  % wrap around circle)
   if ~isempty(zz2)
       p_angos1(zz2) = p_angos1(zz2) + 360;
   end
   
      zz2 = find(p_angos2 > 360);  % wrap around circle)
   if ~isempty(zz2)
       p_angos2(zz2) = p_angos2(zz2) - 360;
   end
   zz2 = find(p_angos2 < 0);  % wrap around circle)
   if ~isempty(zz2)
       p_angos2(zz2) = p_angos2(zz2) + 360;
   end
   
   %*********
   %*** angos are relative to target, put back into absolute coordinates
   angos = angos + tangs;  % offsets from target angle
   zz2 = find(angos > 360);  % wrap around circle)
   if ~isempty(zz2)
       angos(zz2) = angos(zz2) - 360;
   end
   zz2 = find(angos < 0);  % wrap around circle)
   if ~isempty(zz2)
       angos(zz2) = angos(zz2) + 360;
   end
   %*****
   jitter = randn(size(tangs,1),size(tangs,2)) * (0.5 * (180/16));
   atangs = tangs+jitter;
   za = find( atangs > 360);
   atangs(za) = atangs(za)-360;
   za = find( atangs < 0);
   atangs(za) = atangs(za)+360;
   %*************************
   
   %*******   
   if (Targ == 1)
      ha = plot(atangs,angos,'go'); hold on;
      set(ha,'Color',[0,0.5,0]);
      set(ha,'Markersize',2);
      hb = plot(atangs,p_angos,'go'); hold on;
      set(hb,'Color',[1,0.549,0]);
      set(hb,'Markersize',2);
   else
      ha = plot(atangs,angos,'ko'); hold on;
      set(ha,'Color',[0,0,0]);    
      set(ha,'Markersize',2);     
      hb = plot(atangs,[p_angos1; p_angos2],'go'); hold on;
      set(hb,'Color',[1,0.549,0]);
      set(hb,'Markersize',2);
   end
   hb = plot([0,360],[0,360],'k-'); 
   set(hb,'Linewidth',2);
   axis([0 360 0 360]);
   xlabel('True motion direction');
   ylabel('PFR motion direction');
   %*****
   if (Targ == 1)
      title('Saccade Target Motion');
   else
      title('Other Aperture Motion');
   end
   %****
   
return;

function plot_PFR_angles(Info,H)  % plot of angle distribution
          
   subplot(H);
   %**********
   ha = plot(Info.PFR_TargHistX,Info.PFR_TargHist,'g.-'); hold on;
   set(ha,'Color',[0,0.5,0]);
   set(ha,'Linewidth',1);
   hb = plot(Info.PFR_TargHistX,Info.PFR_DistHist,'k.-'); hold on;
   set(hb,'Linewidth',1);
   hc = plot(Info.PFR_TargHistX,Info.PFR_Model_TargHist,'k.-'); hold on;
   set(hc,'Color',[1,0.549,0]);
   set(hc,'Linewidth',1); 
   hd = plot(Info.PFR_TargHistX,Info.PFR_Model_DistHist,'k.-'); hold on;
   set(hd,'Color',[1,0.549,0]);
   set(hd,'Linewidth',1); 
   %***********
   axis tight;
   V = axis;
   axis([-180 180 0 (V(4)*1.2)]);
   xlabel('Angle (degs)');
   ylabel('Probability');
   mederr = Info.PFR_TargErr;
   mederr2 = Info.PFR_DistErr;
   title(sprintf('Median Err( Target: %5.1f, Other: %5.1f)',mederr,mederr2));
   %**********
   
return;


function plot_PFR_vectors(Info,H); % plot of PFR vectors

   subplot(H);
   %****** plot a vector distribution   
   plot(0,0,'ko'); hold on;
   dots = Info.pfrs{1}(:,1); 
   cros = Info.pfrs{1}(:,2); 
   %*** throw out the outliers (huge PFR due to 2nd sac)
   tglength =  Info.OpenLoopSpeed * (Info.OpenLoopEnd-Info.OpenLoopStart);
   mago = abs( dots + i * cros);
   zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
   %***
   for k = 1:length(zz)
     kk = zz(k);
     plot([0,cros(kk)],[0,dots(kk)],'k-');
   end
   udot = mean(dots(zz));
   ucro = mean(cros(zz));
   plot([0,0],[0,tglength],'b-','Linewidth',2);
   plot([0,ucro],[0,udot],'c-','Linewidth',2);
   axis(1.0*[-tglength tglength -tglength tglength]);
   xlabel('x velocity');
   ylabel('y velocity');
   title('PFR vector aligned to target motion');
  
return;

function plot_ProbStat(Info,H); % plot of PFR vectors

  subplot(H);
  colormap;
  imagesc(Info.ProbM2,[0.05 0.8]); hold on;
  colorbar('eastoutside');
  ylabel('Choice');
  xlabel('Last Item ');
  title('Sorted by Freq')
  
return;


function plot_PFR_traces(Info,H); % plot of PFR vectors

  subplot(H);
  uu = nanmean(Info.VelTraces{2,1});
  cu = sum(~isnan(Info.VelTraces{2,1})); % trial counts vary by NaNs
  su = nanstd(Info.VelTraces{2,1}) ./ sqrt(cu);
  %su = nanstd(Info.VelTraces{2,1})/sqrt(size(Info.VelTraces{2,1},1));
  %vtt = 1:length(uu);
  dt = 0.001;  % sampling in secs
  vtt = [(floor(Info.PreStart/dt)-30):-30,1:(floor(Info.PostEnd/dt))];
  h2 = plot(vtt,uu,'r-'); hold on;
  set(h2,'LineWidth',2);
  zn = find(~isnan(uu));
  tta = [vtt(zn) fliplr(vtt(zn))];
  yya = [(uu(zn)+(2*su(zn))) fliplr((uu(zn)-(2*su(zn))))];
  fill(tta,yya,[1,0,0],'FaceAlpha',0.3,'Linestyle','none');
  %*******
  uu = nanmean(Info.VelTraces{2,2});
  cu = sum(~isnan(Info.VelTraces{2,2})); % trial counts vary by NaNs
  su = nanstd(Info.VelTraces{2,2}) ./ sqrt(cu);
  h2 = plot(vtt,uu,'b-'); hold on;
  set(h2,'LineWidth',2);
  zn = find(~isnan(uu));
  tta = [vtt(zn) fliplr(vtt(zn))];
  yya = [(uu(zn)+(2*su(zn))) fliplr((uu(zn)-(2*su(zn))))];
  fill(tta,yya,[0,0,1],'FaceAlpha',0.3,'Linestyle','none');  
  %*****
  plot(vtt,zeros(size(vtt)),'k-');
  axis tight;
  V = axis;
  VA = V(3);
  VB = V(4);
  %****** box for saccade
  aa = [-30,-30,0,0];
  bb = [VA,VB,VB,VA];
  fill(aa,bb,[0.8,0.8,0.7],'FaceAlpha',1.0,'Linestyle','none');
  %***** box for open loop
  OA = floor( Info.OpenLoopStart * 1000);
  OB = floor( Info.OpenLoopEnd * 1000);
  aa = [OA,OA,OB,OB];
  bb = [VA,VB,VB,VA];
  fill(aa,bb,[0.6,0.6,0.6],'FaceAlpha',0.3,'Linestyle','none');
  %********
  xlabel('Time (ms)');
  ylabel('Velocity Gain');
  
return;

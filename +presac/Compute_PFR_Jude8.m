function Info = Compute_PFR_Jude6(Info,Exp,H,plotsteps,trialtype,NormDrift)
%******* 
%******  function Info = Compute_PFR_Jude4(Info,Exp,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          H - if integer, produce plots of the results, if figure
%***                  then gives a cell struct with panel handles to plot
%***          plotsteps - if 1, show trial by trial plots of eye trace and
%***                       of the computed PFR
%***          trialtype - 1, flashing dots with no motion 
%***                      2, moving dots
%***                      3, moving Gabors
%***          NormDrift   - if desired to subtract out non-motion drift
%***
%*** Outputs: Info - it updates fields in the Info struct to include 
%***                  information about trial inclusion parameters
%***

%****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0) 
       %****** setup a multi-panel plot (three panels) 
       hh = figure('position',[100 100 1200 400]); 
       HA = sub('position',[0.10 0.15 0.22 0.70]);
       HB = subplot('position',[0.75 0.15 0.22 0.70]);
       HC = subplot('position',[0.40 0.15 0.25 0.75]);
       plot_PFR_unity(Info,HA,1);   % plot true vs estimate direction, targ
       plot_PFR_unity(Info,HB,2);   % plot true vs estimate direction, other
       plot_PFR_angles(Info,HC,0);  % plot of angle distribution
      
       hh4 = figure('position',[400 200 500 500]); 
       HD = subplot('position',[0.15 0.15 0.70 0.70]);
       plot_PFR_angles(Info,HD,1);  % plot the velocity traces
  
       cohs = unique(Info.PFR_CohList);
       if (length(cohs) >= 2)
           hh5 = figure('position',[600 200 800 400]); 
           HD1 = subplot('position',[0.15 0.15 0.35 0.70]);
           plot_PFR_Coh_angles(Info,HD1,1);  % plot the velocity traces
           HD2 = subplot('position',[0.60 0.15 0.35 0.70]);
           plot_PFR_Coh_curves(Info,HD2,1);  % plot the velocity traces
       end
       
       hh6 = figure('position',[500 200 800 400]);
       HF1 = subplot('position',[0.10 0.15 0.35 0.70]);
       plot_PFR_Gains(Info,HF1,1,1);
       HF2 = subplot('position',[0.60 0.15 0.35 0.70]);
       plot_PFR_Gains(Info,HF2,2,1);
      
       %********************
       hh7 = figure('position',[600 300 500 500]); 
       HG = subplot('position',[0.15 0.15 0.70 0.70]);
       plot_SPFR_traces(Info,HG);  % plot the velocity traces
       
       %********************
       hh2 = figure('position',[300 300 500 500]); 
       HD = subplot('position',[0.15 0.15 0.70 0.70]);
       plot_PFR_traces(Info,HD);  % plot the velocity traces
       
       %********************
       hh3 = figure('position',[300 100 800 800]);
       HE = subplot('position',[0.15 0.15 0.70 0.70]);
       plot_PFR_driftbias(Info,HE);
      
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
          %**********
          z = getframe(hh3);  % current figure
          uname = [Info.pathplot,filesep,'PFR3_',Info.tagname,'.png'];
          disp(sprintf('Storing image of PFR3 graph at %s',uname)); 
          imwrite(z.cdata,uname); % store image for later review
          close(hh3);
       end
    else 
       if ~isempty(H)  % H should be a sub-plot handle otherwise
          plot_PFR_unity(Info,H{1},1);   % plot true vs estimate direction
          plot_PFR_unity(Info,H{2},2);   % plot true vs estimate direction
          plot_PFR_angles(Info,H{3},0);  % plot of angle distribution
          plot_PFR_traces(Info,H{4});  % plot the velocity traces
          plot_PFR_driftbias(Info,H{5});  % plot bias 
       end
    end
    return;
  end

  % define parameters for doing the PFR analysis
  VELTHRESH = 40;  % threshold to flag second saccades
  CURVETHRESH = 1.2;  % threshold to flag curved saccades
  %********
  OpenLoopStart = 0.02;
  OpenLoopEnd = 0.06;  % fast enough to not see other curves diverge
  PostStart = 0.00;
  PostEnd = 0.20;
  PreStart = -0.10;
  PreEnd = 0.00;
  SPostEnd = 0.25;
  SPreStart = -0.05;
  StimSpeed = 6.75;  % would be better to read from Exp file
  Max_Eye_Speed = 20; % if great, NaN from smooth traces
  Max_Eye_Gap = 10;  % NaN +/- this many samples from over Max_Eye
  GSIG = 5;  % smooth the velocity profiles extra
  %*****
  Info.pfr_rt = [];
  NT = Info.NumTargs;
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
  Info.SPreStart = SPreStart;
  Info.SPostEnd = SPostEnd;
  Info.PreStart = PreStart;
  Info.PreEnd = PreEnd;
  %***********
  % plotsteps = 0;  % plot trial by trial to confirm it is working
  %***********
  
  %***** NORMED:  assume each target location causes a default drift
  %************:  run through all trials, store PFR by location and motion
  %************:  compute the weighted average drift at each location
  MT = NT*2;  % two spatial confits
  Info.LocationBias = cell(1,MT);  % each cell a list of x,y vec PFR
  Info.MeanLocBias = cell(1,MT);
  Info.LocationTarg = cell(1,MT);
  Info.SaccadeEnd = cell(1,MT);
  Info.TrialFlag = [];   % perform trial exclusion here in computing net drifts
  Info.SkipTrial = [];
  for k = 1:length( Info.Tlist )
    if ~Info.EPhysList(k) 
        continue;   % skip trial if any problem flagging saccade
    end
    disp(sprintf('Processing Location Bias: trial %d of %d',k,length(Info.Tlist)));
    tr = Info.Tlist(k);  % integer index of an included trial
    
    %**** check if the current trial is the trialtype to include, if not
    %****  continue without including
    Info.SkipTrial(k) = 0;
    if ~isnan(trialtype) 
      if (trialtype == 1)  % flashing dots but no motion
        if ( (Exp.D{tr}.P.motionStimulus ~= 1) || (Exp.D{tr}.P.dotSpeed > 1) )
            Info.SkipTrial(k) = 1;
            continue;
        end
      end
      if (trialtype == 2)  % flashing dots but no motion
        if ( (Exp.D{tr}.P.motionStimulus ~= 1) || (Exp.D{tr}.P.dotSpeed <= 1) ) ...
            || ( ~isfield(Exp.D{tr}.P,'dotLife') || isinf(Exp.D{tr}.P.dotLife)) 
               Info.SkipTrial(k) = 1;
               continue;
        end
      end
      if (trialtype == 4)  % flashing dots but no motion
        if ( (Exp.D{tr}.P.motionStimulus ~= 1) || (Exp.D{tr}.P.dotSpeed <= 1) ) ...
            || ( ~isfield(Exp.D{tr}.P,'dotLife') || ~isinf(Exp.D{tr}.P.dotLife)) 
               Info.SkipTrial(k) = 1;
               continue;
        end
      end
      if (trialtype == 3)  % Gabor stimuli but no motion  
        if ( (Exp.D{tr}.P.motionStimulus ~= 0) || (Exp.D{tr}.P.driftSpeed <= 1) )
            Info.SkipTrial(k) = 1;
            continue;
        end
      end
      %****** 
    end
    %**** given trial, grab eye position data for saccade to aperture    
    eyeSmo = Exp.D{tr}.eyeSmo;  %smoothed eye data (columns: time, x, y)
    Tonsac = Info.SacOnList(k);    % onset time from Info list
    Toffsac = Info.SacOffList(k);  % offset time from Info list
    Tonstim = Info.StimOnList(k);    % onset time from Info list
    if (isnan(Tonsac) || isnan(Toffsac))
        continue;   % skip trial if any problem flagging saccade
    end
    
    
    %***********************************************************
    %********** THIS SUBSECTION IS FOR TRIAL EXCLUSION
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
    %******************
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
      if (length(zstart)>2) && (length(zend)>2)
        for tk = zstart(2):zend(1)
          dist = norm([(xx(tk)-xx(tk-1)),(yy(tk)-yy(tk-1))]);
          pathsum = pathsum + dist;
        end
        curveportion = (pathsum/sacdist);
      else
        curveportion = 1;
      end
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
    if (trialflag ~= 1)
        continue;   % don't include errant trials in drift average
    end
    %************************
    
    %compute PFR vector
    openloop = find((eyeSmo(:,1) >= (Toffsac+OpenLoopStart)) & (eyeSmo(:,1) <= (Toffsac+OpenLoopEnd)));
    ott = eyeSmo(openloop,1);
    oxx = eyeSmo(openloop,2);
    oyy = eyeSmo(openloop,3);
    %compute startposition and open loop vector (openvec)
    if isempty(openloop)
        continue;
    end
    
    %****** plot a regression line (change PFR to that?)
    XX = [ones(size(ott)) ott];
    bx = regress(oxx,XX);
    by = regress(oyy,XX);
    xhat = XX * bx;
    yhat = XX * by;
    newpfr = [xhat(end)-xhat(1),yhat(end)-yhat(1)];  % new PFR
    %***********
    
    startxx = eyeSmo(openloop(1),2);
    startyy = eyeSmo(openloop(1),3);
    startposition = [startxx, startyy];
    % pxx = (eyeSmo(openloop(end),2) - startxx); % x position difference
    % pyy = (eyeSmo(openloop(end),3) - startyy); % y position difference
    pxx = newpfr(1);
    pyy = newpfr(2);
    
    %******* store vectors, is there a net drift per location
    tg = Info.TargList(k,1); % which of 3 apertures was the target 
    if (isfield(Exp.D{tr}.PR,'SpatialConfig'))
      spcon = Exp.D{tr}.PR.SpatialConfig;
    else
      if (0) % isfield(Exp.D{tr}.P,'spatialconfig')
        spcon = Exp.D{tr}.P.spatialconfig;
      else
        if (Exp.D{tr}.P.RF_Y < 0)  
         spcon = 1;
        else
         spcon = 2;
        end
      end
    end
    itg = (spcon-1)*NT + tg;
    if isempty(Info.LocationTarg{itg})
       posx = Info.TargPosX(k,tg);
       posy = Info.TargPosY(k,tg);
       Info.LocationTarg{itg} = [posx,posy];
    end
    Info.SaccadeEnd{itg} = [ Info.SaccadeEnd{itg} ; [xend,yend]];
    %********
    if isnan(pxx) || isnan(pyy)    % eye tracker failed
        Info.TrialFlag(k) = 5;   % failure in eye track
        continue;
    end
    Info.LocationBias{itg} = [Info.LocationBias{itg} ; [pxx,pyy] ];
 end
 for k = 1:MT
      if isempty(Info.LocationBias{k})
         Info.MeanLocBias{k} = [0,0];
      else 
         if length(Info.LocationBias{k}) <= 2
           Info.MeanLocBias{k} = [0,0];     
         else
           Info.MeanLocBias{k} = nanmean(Info.LocationBias{k});  % use median, reduce outlier?
         end
      end
 end
 %*************  
 
 %***** Compute corrections for camera rotation
 angdev = 0;
 angcnt = 0;
 for k = 1:size(Info.LocationTarg,2)
     if ~isempty(Info.LocationTarg{k})
         %******* compute dots and cross-prods per end-point
         txx = Info.LocationTarg{k}(1);
         tyy = Info.LocationTarg{k}(2);
         nomo = norm([txx,tyy]);
         for kk = 1:size(Info.SaccadeEnd{k},1)
            sxx = Info.SaccadeEnd{k}(kk,1);
            syy = Info.SaccadeEnd{k}(kk,2);
            somo = norm([sxx,syy]);
            %*********
            doto = ((txx*sxx)+(tyy*syy))/(somo*nomo);
            cros = ((txx*syy)-(tyy*sxx))/(somo*nomo);
            ango = angle( complex(doto,cros));
            %*******
            angcnt = angcnt + 1;
            angdev = angdev + ango;
         end
     end
 end
 mango = mean(ango);
 mdeg = mango*(180/pi);
 disp(sprintf('Correcting median camera rotation of %5.1f degs',mdeg));
 %*********
 Info.Mroto = [[cos(mango) sin(mango)]; [-sin(mango) cos(mango)]];
 %***************
 
%  figure(10);
%  for k = 1:size(Info.LocationTarg,2)
%      if ~isempty(Info.LocationTarg{k})
%          
%          plot(Info.SaccadeEnd{k}(:,1)',Info.SaccadeEnd{k}(:,2)','g.'); hold on;
%          
%          sacrot = (Info.Mroto * Info.SaccadeEnd{k}')';
%          
%          plot(sacrot(:,1)',sacrot(:,2)','m.'); hold on;
%          
%          %******* plot correction **********
%          
%          if (k <= 3)
%            plot(Info.LocationTarg{k}(1),Info.LocationTarg{k}(2),'k+');
%          else
%            plot(Info.LocationTarg{k}(1),Info.LocationTarg{k}(2),'b+');             
%          end
%      end
%      axis([-8 8 -8 8]);
%  end
%  input('check points');
%  
 
  %******** second pass through data will normalize for mean drift
  if (plotsteps == 1)
            %**********plot some of the data
            h = figure(1); 
            set(h,'position',[100 100 1200 800]);  hold off;
  end
  
  %***** store velocity traces
  Info.VelTraces = cell(2,3);
  Info.SVelTraces = cell(1,3);
  % loop through the list of trials accepted for analyses
  disp('Computing PFR statistics ...');
  Info.PFR_DotList = nan(size(Info.Tlist));  % store dot prods 
  Info.PFR_CohList = [];
  Info.PFR_MotiList = [];
  Info.PFR_TangList = [];
  Info.PFR_GainList = [];
  Info.PFR_FomoList = [];
  Info.PFR_VelTangList = cell(2,3);  % to match vel traces stored
  for k = 1:length( Info.Tlist )
    if ~Info.EPhysList(k) 
        continue;   % skip trial if any problem flagging saccade
    end
    if Info.SkipTrial(k)  % if not the right trialtype, skip
        continue;
    end
    %*********
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
    %**************************
     
    %compute PFR vector
    openloop = find((eyeSmo(:,1) >= (Toffsac+OpenLoopStart)) & (eyeSmo(:,1) <= (Toffsac+OpenLoopEnd)));
    ott = eyeSmo(openloop,1);
    oxx = eyeSmo(openloop,2);
    oyy = eyeSmo(openloop,3);
   
    %compute startposition and open loop vector (openvec)
    if isempty(openloop)
        continue;
    end
   
    %*********
    XX = [ones(size(ott)) ott];
    bx = regress(oxx,XX);
    by = regress(oyy,XX);
    xhat = XX * bx;
    yhat = XX * by;
    newpfr = [xhat(end)-xhat(1),yhat(end)-yhat(1)];  % new PFR
    %***********
    
    startxx = eyeSmo(openloop(1),2);
    startyy = eyeSmo(openloop(1),3);
    startposition = [startxx, startyy];    
    % pxx = (eyeSmo(openloop(end),2) - startxx); % x position difference
    % pyy = (eyeSmo(openloop(end),3) - startyy); % y position difference
    pxx = newpfr(1);
    pyy = newpfr(2);
    % openvec = [pxx,pyy];
   
    %******* determine mean drift per target location
    tg = Info.TargList(k,1); % which of 3 apertures was the target
    if (isfield(Exp.D{tr}.PR,'SpatialConfig'))
      spcon = Exp.D{tr}.PR.SpatialConfig;
    else
      if (0) % isfield(Exp.D{tr}.P,'spatialconfig')
        spcon = Exp.D{tr}.P.spatialconfig;
      else
        if (Exp.D{tr}.P.RF_Y < 0)  
         spcon = 1;
        else
         spcon = 2;
        end
      end
    end
    itg = (spcon-1)*NT + tg;
    mxx = Info.MeanLocBias{itg}(1);
    myy = Info.MeanLocBias{itg}(2);  % median PFR shift per location
    vmxx = mxx/(OpenLoopEnd-OpenLoopStart);
    vmyy = myy/(OpenLoopEnd-OpenLoopStart); % median vel vector per loc
    %*********
    if (NormDrift)
        openvec = [(pxx-mxx),(pyy-myy)];
    else
        openvec = [pxx,pyy];
    end
    
    %***** correct rotation (obligatory)
    openvec = ( Info.Mroto * openvec' )';
    %********************   
    
    %**** compute tangent angle and store per trial
    tg = Info.TargList(k,1); % which of 3 apertures was the target 
    tpos = [Info.TargPosX(k,tg),Info.TargPosY(k,tg)];
    tpos = tpos / norm(tpos);
    tgmotion = Info.TargMotion(k,tg);  % motion in degrees
    ango = tgmotion * pi/180;
    tgvec = [cos(ango), sin(ango)];
    doto = (tgvec(1) * tpos(1)) + (tgvec(2) * tpos(2));
    croo = -((tgvec(1) * tpos(2)) - (tgvec(2) * tpos(1))); 
    veltang = angle( doto +  i * croo ) * (180/pi);
    %****** Store the angle somewhere
    Info.PFR_TangList = [Info.PFR_TangList ; veltang];  
    Info.PFR_MotiList = [Info.PFR_MotiList ; (ango*180/pi)];
    %***** compute PFR dot prod and store
    dotsp = Exp.D{tr}.P.dotSpeed;  %smoothed eye data (column
    if (StimSpeed ~= dotsp)
        StimSpeed = dotsp;
    end
    endspeed = StimSpeed * (OpenLoopEnd-OpenLoopStart);
    pfrg = ((openvec(1)*tgvec(1))+(openvec(2)*tgvec(2))) / endspeed;
    Info.PFR_GainList = [Info.PFR_GainList ; pfrg];
    if isnan( Exp.D{tr}.PR.changori )  % if NaN, it blanked
        Info.PFR_FomoList = [Info.PFR_FomoList ; 0]; 
    else
       if (Exp.D{tr}.PR.changori == Exp.D{tr}.PR.targori) %here stayed the same
          Info.PFR_FomoList = [Info.PFR_FomoList ; 1]; 
       else
          Info.PFR_FomoList = [Info.PFR_FomoList ; NaN];       
       end
    end
    %****
    
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
  
    %******* compute stim locked vel trace
    ploop3 = find((eyeSmo(:,1) >= (Tonstim+SPreStart)) & (eyeSmo(:,1) <= (Tonstim+SPostEnd)));
    if isempty(ploop3)
        ploop3 = Inf;  % will NaN out ...
    end
    plop3 = ploop3(1) + (0:floor((SPostEnd-SPreStart)/0.001));
    if (max(plop3) > size(eyeSmo,1) )
       plxx3 = NaN(1,length(plop3));
       plyy3 = NaN(1,length(plop3));
       pltt3 = 1:length(plop3);
       disp('Failed to reach limit of eye position');
    else    
       plxx3 = eyeSmo(plop3,2);
       plyy3 = eyeSmo(plop3,3);
       pltt3 = eyeSmo(plop3,1);
    end
    %***** compute the velocity traces
    svxx = [];
    svyy = [];
    svsp = [];
    svtt = [];
    for t = 2:length(plxx3)
       svelx = (plxx3(t) - plxx3(t-1))/(pltt3(t)-pltt3(t-1));  % dva / secs
       svely = (plyy3(t) - plyy3(t-1))/(pltt3(t)-pltt3(t-1));  % dva / secs
       svabs = norm([svelx,svely]);
       svtt = [svtt ; ((pltt3(t)+pltt3(t-1))/2)];  % mid-point of time between the two points used to compute velocity
       svxx = [svxx svelx];
       svyy = [svyy svely];
       svsp = [svsp svabs];
    end
    %******* Apply additional smoothing for velocity
    svxx = spikestats.gauss_smooth(svxx,GSIG);
    svyy = spikestats.gauss_smooth(svyy,GSIG);
    svsp = spikestats.gauss_smooth(svsp,GSIG);
    prelip = 0.02;  % 10 ms min before saccade starts
    if ((Tonsac-prelip) < (Tonstim+SPostEnd))
       dropit = floor( ((Tonstim+SPostEnd)-(Tonsac-prelip))/0.001);
       svxx((end-dropit):end) = NaN;
       svyy((end-dropit):end) = NaN;
       svsp((end-dropit):end) = NaN;
    end
    %****** find and delete any velocity violation trials
    for t = 1:length(svsp)
        if ~isnan(svsp(t)) && (svsp(t) > Max_Eye_Speed)
           ta = max(1,(t-Max_Eye_Gap));
           tb = length(svsp); % min(length(nvsp),(t+Max_Eye_Gap));
           svxx(ta:tb) = NaN; % throw out rest of trial if saccade
           svyy(ta:tb) = NaN;
        end
    end
   
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
    %**** do we need to project away mean drift?
    %*** compute same but relative to other targets 
    ovdot1 = ((tgvec2(1) * nvxx) + (tgvec2(2) * nvyy))/(tglength * StimSpeed);
    ovdot2 = ((tgvec3(1) * nvxx) + (tgvec3(2) * nvyy))/(tglength * StimSpeed);
    %********
    svdot = ((tgvec(1) * svxx) + (tgvec(2) * svyy))/(tglength * StimSpeed);
    %**** do we need to project away mean drift?
    %*** compute same but relative to other targets 
    sovdot1 = ((tgvec2(1) * svxx) + (tgvec2(2) * svyy))/(tglength * StimSpeed);
    sovdot2 = ((tgvec3(1) * svxx) + (tgvec3(2) * svyy))/(tglength * StimSpeed);
    %***********
    if isfield(Exp.D{tr}.PR,'targori')
         if isnan( Exp.D{tr}.PR.changori )  % if NaN, it blanked
               Info.VelTraces{1,2} = [Info.VelTraces{1,2} ; vdot];
               Info.VelTraces{2,2} = [Info.VelTraces{2,2} ; nvdot];
               Info.PFR_VelTangList{1,2} = [Info.PFR_VelTangList{1,2} ; veltang];
               Info.PFR_VelTangList{2,2} = [Info.PFR_VelTangList{2,2} ; veltang];
         else
            if (Exp.D{tr}.PR.changori == Exp.D{tr}.PR.targori) %here stayed the same
               Info.VelTraces{1,1} = [Info.VelTraces{1,1} ; vdot];
               Info.VelTraces{2,1} = [Info.VelTraces{2,1} ; nvdot];
               Info.PFR_VelTangList{1,1} = [Info.PFR_VelTangList{1,1} ; veltang];
               Info.PFR_VelTangList{2,1} = [Info.PFR_VelTangList{2,1} ; veltang];
            end
         end
         %***** relative to other two targets
         Info.VelTraces{1,3} = [Info.VelTraces{1,3} ; ovdot1];
         Info.VelTraces{2,3} = [Info.VelTraces{2,3} ; ovdot2];
         Info.PFR_VelTangList{1,3} = [Info.PFR_VelTangList{1,3} ; 0];  % does not apply for other
         Info.PFR_VelTangList{2,3} = [Info.PFR_VelTangList{2,3} ; 0];
         %*************
         Info.SVelTraces{1,1} = [Info.SVelTraces{1,1} ; svdot];
         Info.SVelTraces{1,2} = [Info.SVelTraces{1,2} ; sovdot1];
         Info.SVelTraces{1,3} = [Info.SVelTraces{1,3} ; sovdot2];
         %*******   
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
    %*****
    if isfield(Info,'CohList')
        Info.PFR_CohList = [Info.PFR_CohList ; Info.CohList(k)];
    else
        Info.PFR_CohList = [];
    end
    
    %****** look at trial flag, and decode for graphing it ******
    %****** the computation of flagged trials is above
    secondsac = 0;
    curvesac = 0;
    if (Info.TrialFlag(k) == 2)
        secondsac = 1;
    end
    if (Info.TrialFlag(k) == 3)
        curvesac = 1;
    end
    
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
            %****** plot a regression line (change PFR to that?)
            %  XX = [ones(size(ott)) ott];
            %  bx = regress(oxx,XX);
            %  by = regress(oyy,XX);
            %  xhat = XX * bx;
            %  yhat = XX * by;
            %  newpfr = [xhat(end)-xhat(1),yhat(end)-yhat(1)];  % new PFR
            %***********
            
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
              %*******
              plot(xhat,yhat,'r-');
              %********************
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
            if (0)
              wid = 2;
              xa = Info.TargPosX(k,tg);
              ya = Info.TargPosY(k,tg);
              axis([xa-wid xa+wid ya-wid ya+wid]);  % bound in -10 to +10 window
              wid = 1.95;
              plot([(xa-wid),(xa-wid),(xa+wid),(xa+wid),(xa-wid)],...
                 [(ya-wid),(ya+wid),(ya+wid),(ya-wid),(ya-wid)],'g-');
              title(sprintf('2ndsac(%d) Curve(%d)',secondsac,curvesac));
            end
            
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
  
  %***** then do a Jacknife for error bars
  dhist = [];
  LA = length(angos);
  JN = 10;
  for jk = 1:JN
       aset = [1:(1+(jk-1)*floor(LA/JN)),((jk+1)*floor(LA/JN)):LA];
       vvy = hist(angos(aset),vx);
       vvy = vvy / sum(vvy);
       dhist = [dhist ; vvy];
  end
  Info.PFR_TargHistSU = std(dhist)*sqrt(JN-1);
 
  if ~isempty(Info.PFR_CohList)
    cohs = unique(Info.PFR_CohList(zz));
    if (length(cohs) > 1)
      Info.PFR_CohHist = cell(2,length(cohs));
      Info.PFR_MedCohHist = NaN(2,length(cohs));
      Info.PFR_StdCohHist = NaN(2,length(cohs));
      Info.PFR_Shaun_ProCor = NaN(2,length(cohs));
      %*** build vectors for Shaun analysis too
      shaunY = []; % PFR per trial
      shaunS = []; % actual direction on trial
      shaunBW = []; % integer coherence per trial
      %*******
      for ck = 1:length(cohs)
       cobo = cohs(ck);
       za = find( Info.PFR_CohList(zz) == cobo);
       zangos = angle( dots(zz(za)) +  i * cros(zz(za)) ) * (180/pi);  % go from radians to degs
       zvy = hist(zangos,vx);
       zvy = zvy / sum(zvy);   % normalize
       Info.PFR_CohHist{1,ck} = zvy;
       %***
       zza = find( abs(zangos) <= 27);
       Info.PFR_Shaun_ProCor(1,ck) = length(zza)/length(zangos);
       %***** build for Shaun
       shaunY = [shaunY ; zangos];  % errors are from -180 to + 180
       shaunS = [shaunS ; zeros(size(zangos))];
       shaunBW = [shaunBW ; (ck*ones(size(zangos)))];
       %*******
       Info.PFR_MedCohHist(1,ck) = median(abs(zangos));
       Info.PFR_StdCohHist(1,ck) = std(zangos);
       %********* Jacknife
       dhist = [];
       dmed = [];
       dstd = [];
       pcor = [];
       LA = length(zangos);
       for jk = 1:JN
          aset = [1:(1+(jk-1)*floor(LA/JN)),((jk+1)*floor(LA/JN)):LA];
          vvy = hist(zangos(aset),vx);
          vvy = vvy / sum(vvy);
          dhist = [dhist ; vvy];
          %******
          zza = find( abs(zangos(aset)) <= 27);
          pcor = [pcor ; (length(zza)/length(zangos(aset)))];
          %*********
          dmed = [dmed ; median(abs(zangos(aset)))];
          dstd = [dstd ; std(zangos(aset))];
          %*********
       end
       Info.PFR_CohHist{2,ck} = std(dhist)*sqrt(JN-1);
       Info.PFR_MedCohHist(2,ck) = std(dmed)*sqrt(JN-1);
       Info.PFR_StdCohHist(2,ck) = std(dstd)*sqrt(JN-1);
       Info.PFR_Shaun_ProCor(2,ck) = std(pcor)*sqrt(JN-1);
       %*******************
      end
      %******** Use Shaun's model fitting
      Info.PFR_Shaun_StdCohHist = NaN(2,length(cohs));
      Info.PFR_Shaun_Lapse = NaN(2,1);
      %***** 
      nBw = length(cohs);
      m = wn_bds_ulapse(nBw);
      p = m.fit(shaunY,shaunS,shaunBW)';
      Info.PFR_Shaun_StdCohHist(1,:) = [p(2:end)]; %p(1) is lapse rate
      Info.PFR_Shaun_Lapse(1,1) = p(1);
      %*** get error bars via Jacknife
      jackp = [];
      LA = length(shaunY);
      rn = randperm(LA);
      for jk = 1:JN
          aset = [1:(1+(jk-1)*floor(LA/JN)),((jk+1)*floor(LA/JN)):LA];
          raset = rn(aset);
          jp = m.fit(shaunY(raset),shaunS(raset),shaunBW(raset));
          jackp = [jackp ; jp'];
       end
       p = std(jackp)*sqrt(JN-1);
       Info.PFR_Shaun_StdCohHist(2,:) = [p(2:end)]; %p(1) is lapse rate
       Info.PFR_Shaun_Lapse(2,1) = p(1);
       %********    
       Info.PFR_Shaun_StdCohHist
       Info.PFR_Shaun_ProCor
    end
  end
  
  %****** Gains ****
  [xtang,uu,su] = plot_PFR_Gains(Info,[],1,0);
  Info.PFR_TangGain = [xtang' uu su];
  [xtang,uu,su] = plot_PFR_Gains(Info,[],2,0);
  Info.PFR_MotiGain = [xtang' uu su];
  
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
  %***** then do a Jacknife for error bars
  dhist = [];
  LA = length(angos);
  JN = 10;
  for jk = 1:JN
       aset = [1:(1+(jk-1)*floor(LA/JN)),((jk+1)*floor(LA/JN)):LA];
       vvy = hist(angos(aset),vx);
       vvy = vvy / sum(vvy);
       dhist = [dhist ; vvy];
  end
  Info.PFR_DistHistSU = std(dhist)*sqrt(JN-1);
   
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
  Info.VelTraces_U3 = 0.5 * ( nanmean(Info.VelTraces{1,3}) + nanmean(Info.VelTraces{2,3}));
  Info.VelTraces_S3 = 0.5 * ( nanstd(Info.VelTraces{1,3}) + nanstd(Info.VelTraces{2,3}) ) / ...
                             sqrt( size(Info.VelTraces{1,3},1) + size(Info.VelTraces{2,3},1));
  %*******
  Info.SVelTraces_U1= nanmean(Info.SVelTraces{1,1});
  Info.SVelTraces_S1 = nanstd(Info.SVelTraces{1,1})/sqrt(size(Info.SVelTraces{1,1},1));
  Info.SVelTraces_U3 = 0.5 * ( nanmean(Info.SVelTraces{1,2}) + nanmean(Info.SVelTraces{1,3}));
  Info.SVelTraces_S3 = 0.5 * ( nanstd(Info.SVelTraces{1,2}) + nanstd(Info.SVelTraces{1,3}) ) / ...
                                sqrt( size(Info.SVelTraces{1,2},1) + size(Info.SVelTraces{1,3},1));
                         
  %******* Same, but tangent only  *************
  Info.VelTangThresh = 30;
  Info.VelTangThreshA = 90-Info.VelTangThresh;  % tangent trials from sac
  Info.VelTangThreshB = 90+Info.VelTangThresh;
  zn21 = find( (abs(Info.PFR_VelTangList{2,1}) >= Info.VelTangThreshA) & ...
               (abs(Info.PFR_VelTangList{2,1}) <= Info.VelTangThreshB) );
  zn22 = find( (abs(Info.PFR_VelTangList{2,2}) >= Info.VelTangThreshA) & ...
               (abs(Info.PFR_VelTangList{2,2}) <= Info.VelTangThreshB) );
  Info.VelTraces_ZU1= nanmean(Info.VelTraces{2,1}(zn21,:));
  Info.VelTraces_ZS1 = nanstd(Info.VelTraces{2,1}(zn21,:))/sqrt(size(Info.VelTraces{2,1}(zn21,:),1));
  Info.VelTraces_ZU2 = nanmean(Info.VelTraces{2,2}(zn22,:));
  Info.VelTraces_ZS2 = nanstd(Info.VelTraces{2,2}(zn22,:))/sqrt(size(Info.VelTraces{2,2}(zn22,:),1));                       
  %**************
  
  %******** last step is to plot results if desired
  disp(sprintf('PFR statistics finished computing: Targ(%5.2f) Ign(%5.2f)',...
          Info.PFR_TargErr,Info.PFR_DistErr));

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
   angos = angle( dots +  i * cros ) *  (180/pi); %(180/pi);  % go from radians to degs
   tangs = tangs; % match vectors
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
   %*******   
   if (Targ == 1)
      ha = plot(atangs,angos,'go'); hold on;
      set(ha,'Color',[0,0.5,0]);
      set(ha,'Markersize',2);
   else
      ha = plot(atangs,angos,'ko'); hold on;
      set(ha,'Color',[0,0,0]);    
      set(ha,'Markersize',2);
   end
   hb = plot([0,360],[0,360],'k-'); 
   set(hb,'Linewidth',2);
   axis([0 360 0 360]);
   xlabel('True motion direction');
   ylabel('PFR motion direction');
   %*****
   if (Targ == 1)
      title('Saccade Targ Motion');
   else
      title('Other Aperture Motion');
   end
   %****
   
return;

function plot_PFR_angles(Info,H,CrossValid)  % plot of angle distribution
          
   subplot(H);
   %**********
   ha = plot(Info.PFR_TargHistX,Info.PFR_TargHist,'g.-'); hold on;
   set(ha,'Color',[0,0.5,0]);
   if (CrossValid)
       tta = [Info.PFR_TargHistX fliplr(Info.PFR_TargHistX)];
       uu = Info.PFR_TargHist;
       su = Info.PFR_TargHistSU;
       yya = [(uu+(2*su)) fliplr((uu-(2*su)))];
       fill(tta,yya,[0,0.5,0],'FaceAlpha',0.3,'Linestyle','none');
   end
   set(ha,'Linewidth',2);
   hb = plot(Info.PFR_TargHistX,Info.PFR_DistHist,'k.-'); hold on;
   set(hb,'Linewidth',2);
   if (CrossValid)
       tta = [Info.PFR_TargHistX fliplr(Info.PFR_TargHistX)];
       uu = Info.PFR_DistHist;
       su = Info.PFR_DistHistSU;
       yya = [(uu+(2*su)) fliplr((uu-(2*su)))];
       fill(tta,yya,[0,0,0],'FaceAlpha',0.3,'Linestyle','none');      
   end
   %***********
   axis tight;
   V = axis;
   axis([-180 180 0 0.20]);
   plot([0,0],[0,0.20],'k-');
   xlabel('Angle (degs)');
   ylabel('Probability');
   mederr = Info.PFR_TargErr;
   mederr2 = Info.PFR_DistErr;
   title(sprintf('Median Err( Targ: %5.1f, Other: %5.1f)',mederr,mederr2));
   %**********
   
return;


function plot_PFR_Coh_angles(Info,H,CrossValid)  % plot of angle distribution          
   subplot(H);
   %**********
   cohs = unique(Info.PFR_CohList);
   for ck = 1:length(cohs)
     val = ((ck-1)/length(cohs));
     colo = [0,(1-val),0];
     %*****
     ha = plot(Info.PFR_TargHistX,Info.PFR_CohHist{1,ck},'g-'); hold on;
     set(ha,'Color',colo);
     set(ha,'LineWidth',2);
     if (CrossValid)
       tta = [Info.PFR_TargHistX fliplr(Info.PFR_TargHistX)];
       uu = Info.PFR_CohHist{1,ck};
       su = Info.PFR_CohHist{2,ck};
       yya = [(uu+(2*su)) fliplr((uu-(2*su)))];
       fill(tta,yya,colo,'FaceAlpha',0.3,'Linestyle','none');
     end
   end
   %***********
   axis tight;
   V = axis;
   axis([-180 180 0 0.20]);
   plot([0,0],[0,0.20],'k-');
   xlabel('Angle (degs)');
   ylabel('Probability');
   title(sprintf('By Coherence (Green highest)'));
   %**********   
return;

function plot_PFR_Coh_curves(Info,H,CrossValid)  % plot of angle distribution          
   subplot(H);
   %**********
   colo = [0,0.5,0];
   icolo = [0,1,0];
   cohs = unique(Info.PFR_CohList);
   %*******
   uu = Info.PFR_StdCohHist(1,:);
   su = Info.PFR_StdCohHist(2,:);
   tta = ((360 - cohs)/360)';   
   ha = plot(tta,uu,'g.-'); hold on;
   set(ha,'Color',colo);
   set(ha,'LineWidth',2);
   set(ha,'Markersize',20);
   if (CrossValid)
     yya = [(uu+(2*su)) fliplr((uu-(2*su)))];
     ttb = [tta fliplr(tta)];
     fill(ttb,yya,icolo,'FaceAlpha',0.3,'Linestyle','none');
   end
   %*******
   colo = [0,0,0.5];
   icolo = [0,0,1];
   uu = Info.PFR_MedCohHist(1,:);
   su = Info.PFR_MedCohHist(2,:);
   tta = ((360 - cohs)/360)';   
   ha = plot(tta,uu,'k.-'); hold on;
   set(ha,'Color',colo);
   set(ha,'LineWidth',2);
   set(ha,'Markersize',20);
   if (CrossValid)
     yya = [(uu+(2*su)) fliplr((uu-(2*su)))];
     ttb = [tta fliplr(tta)];
     fill(ttb,yya,icolo,'FaceAlpha',0.3,'Linestyle','none');
   end
   %*******
   colo = [0.5,0,0];
   icolo = [1,0,0];
   uu = Info.PFR_Shaun_StdCohHist(1,:);
   su = Info.PFR_Shaun_StdCohHist(2,:);
   tta = ((360 - cohs)/360)';   
   ha = plot(tta,uu,'k.-'); hold on;
   set(ha,'Color',colo);
   set(ha,'LineWidth',2);
   set(ha,'Markersize',20);
   if (CrossValid)
     yya = [(uu+(2*su)) fliplr((uu-(2*su)))];
     ttb = [tta fliplr(tta)];
     fill(ttb,yya,icolo,'FaceAlpha',0.3,'Linestyle','none');
   end
   %***********
   set(gca,'xtick',0:0.25:1);
   set(gca,'ytick',0:45:180);
   axis([0 1 0 180]);
   box off;
   xlabel('Signal Strength (SS)');
   ylabel('Angular Error (Std)');
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
  %***********
  uu = 0.5 * ( nanmean(Info.VelTraces{1,3}) + nanmean(Info.VelTraces{2,3}));
  cu = sum(~isnan(Info.VelTraces{1,3})) + sum(~isnan(Info.VelTraces{2,3})); 
  su = 0.5 * ( nanstd(Info.VelTraces{1,3}) + nanstd(Info.VelTraces{2,3}));
  su = su ./ sqrt( cu );
  h2 = plot(vtt,uu,'k-'); hold on;
  set(h2,'LineWidth',2);
  zn = find(~isnan(uu));
  tta = [vtt(zn) fliplr(vtt(zn))];
  yya = [(uu(zn)+(2*su(zn))) fliplr((uu(zn)-(2*su(zn))))];
  fill(tta,yya,[0,0,0],'FaceAlpha',0.3,'Linestyle','none');  
  %*****
  plot(vtt,zeros(size(vtt)),'k-');
  axis tight;
  V = axis;
  VA = V(3);
  VB = V(4);
  VA = -0.05;
  VB = 0.45;
  axis([V(1) V(2) VA VB]);
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

function plot_SPFR_traces(Info,H)  % plot of PFR vectors

  subplot(H);
  uu = nanmean(Info.SVelTraces{1,1});
  cu = sum(~isnan(Info.SVelTraces{1,1})); % trial counts vary by NaNs
  su = nanstd(Info.SVelTraces{1,1}) ./ sqrt(cu);
  dt = 0.001;  % sampling in secs
  vtt = [(floor(Info.SPreStart/dt)+1):(floor(Info.SPostEnd/dt))];
  h2 = plot(vtt,uu,'r-'); hold on;
  set(h2,'LineWidth',2);
  zn = find(~isnan(uu));
  tta = [vtt(zn) fliplr(vtt(zn))];
  yya = [(uu(zn)+(2*su(zn))) fliplr((uu(zn)-(2*su(zn))))];
  fill(tta,yya,[1,0,0],'FaceAlpha',0.3,'Linestyle','none');
  %***********
  uu = 0.5 * ( nanmean(Info.SVelTraces{1,2}) + nanmean(Info.SVelTraces{1,3}));
  cu = sum(~isnan(Info.SVelTraces{1,2})) + sum(~isnan(Info.SVelTraces{1,3})); 
  su = 0.5 * ( nanstd(Info.SVelTraces{1,2}) + nanstd(Info.SVelTraces{1,3}));
  su = su ./ sqrt( cu );
  h2 = plot(vtt,uu,'k-'); hold on;
  set(h2,'LineWidth',2);
  zn = find(~isnan(uu));
  tta = [vtt(zn) fliplr(vtt(zn))];
  yya = [(uu(zn)+(2*su(zn))) fliplr((uu(zn)-(2*su(zn))))];
  fill(tta,yya,[0,0,0],'FaceAlpha',0.3,'Linestyle','none');  
  %*****
  plot(vtt,zeros(size(vtt)),'k-');
  axis tight;
  V = axis;
  VA = V(3);
  VB = V(4);
  VA = -0.05;
  VB = 0.45;
  axis([V(1) V(2) VA VB]);
  %********
  xlabel('Time (ms)');
  ylabel('Velocity Gain');
  
return;


function [xtang,uu,su] = plot_PFR_Gains(Info,HF1,type,plotit)
  
  if (plotit)
    subplot(HF1);
  end
  gains = Info.PFR_GainList;
  if (type == 1)
     tangos = Info.PFR_TangList;
  else
     tangos = Info.PFR_MotiList;    
  end
  utang = unique(tangos);
  gap = (360/16);
  gapt = (gap*1);
  uu = [];
  su = [];
  xtang = 0:gap:(360-gap);
  for ztang = xtang
     itang = tangos - ztang;
     zz = find( itang > 180);
     itang(zz) = itang(zz) - 360;
     zz = find( itang < -180);
     itang(zz) = itang(zz) + 360;
     zn = find( abs(itang) <= gapt);
     %****
     uval = nanmean( gains(zn) );
     sval = nanstd( gains(zn) )/sqrt(length(zn));
     uu = [uu ; uval];
     su = [su ; sval];
  end
  if (type == 1)
      xtang = xtang + 90;  % reference to vertical
  end
  if (plotit)
    mypolar(xtang,uu,su,[0,0,1]);
    if (type == 1)
      title('Gain by Tangent'); 
    else
      title('Gain by Motion'); 
    end
  end
  
return


function plot_PFR_driftbias(Info,HE)

%******* plot individual drifts per locations in space available
subplot(HE);
NT = length(Info.LocationBias);
MT = floor(NT/2);
%******* find a way to scale biases to fit better for viewing

%********
varo = 0;
for k = 1:NT
    dists = sqrt( sum( (Info.LocationBias{k}') .^ 2) );
    varo = varo + mean(dists);
end
Base = 2;
scalo = Base/varo;
dv = (Base*varo);
dv2 = dv/scalo;
%**********
for k = 1:NT
  if isempty(Info.LocationTarg{k})
      continue;
  end
  posx = Info.LocationTarg{k}(1);
  posy = Info.LocationTarg{k}(2);
  for j = 1:size(Info.LocationBias{k},1)
      pxx = scalo * Info.LocationBias{k}(j,1);
      pyy = scalo * Info.LocationBias{k}(j,2);
      h = plot([posx+0,posx+pxx],[posy+0,posy+pyy],'k-'); hold on;
      set(h,'Color',[0.5,0.5,0.5]);
  end
  mxx = scalo * Info.MeanLocBias{k}(1);
  myy = scalo * Info.MeanLocBias{k}(2);
  h = plot([posx+0,posx+mxx],[posy+0,posy+myy],'b-');
  set(h,'Linewidth',2);
  h = plot(posx,posy,'k.');
  set(h,'Markersize',15);
  %****** draw a box around position
  hc = plot([(posx-dv),(posx+dv),(posx+dv),(posx-dv),(posx-dv)],...
       [(posy-dv),(posy-dv),(posy+dv),(posy+dv),(posy-dv)],'c-'); hold on;
  if (k <= MT)
    set(hc,'Color',[0,0.6,0.6]);
  else
    set(hc,'Color',[0.6,0.6,0.0]);      
  end
  hc = plot([(posx-dv2),(posx+dv2),(posx+dv2),(posx-dv2),(posx-dv2)],...
       [(posy-dv2),(posy-dv2),(posy+dv2),(posy+dv2),(posy-dv2)],'c--'); hold on;
  if (k <= MT)
     set(hc,'Color',[0,0.6,0.6]);
  else
     set(hc,'Color',[0.6,0.6,0.0]);     
  end
  %********
end
axis tight;
V = axis;
maxo = max(abs(V))*1.1;
axis([-maxo maxo -maxo maxo]);
hb = plot(0,0,'k+');
set(hb,'Markersize',15);
xlabel('X position (degs)');
ylabel('Y position (degs)');
title(sprintf('%s : Drifts (Actual Box Side %3.1f)',Info.tagname,((dv*2)/scalo)));
%*********

return;


function mypolar(xtang,uu,su,colo)
    uux = [];
    uuy = [];
    huux = [];
    huuy = [];
    luux = [];
    luuy = [];
    for k = 1:(length(xtang)+1) 
       if (k > length(xtang))
           kk = 1;
       else
           kk = k;
       end
       ango = xtang(kk)*(pi/180);
       uux = [uux ; (cos(ango)*uu(kk))];
       uuy = [uuy ; (sin(ango)*uu(kk))];
       huux = [huux ; (cos(ango)*(uu(kk)+(2*su(kk))))];
       huuy = [huuy ; (sin(ango)*(uu(kk)+(2*su(kk))))];
       luux = [luux ; (cos(ango)*(uu(kk)-(2*su(kk))))];
       luuy = [luuy ; (sin(ango)*(uu(kk)-(2*su(kk))))];
    end
    %**********
    h1 = plot(uux,uuy,'ko-'); hold on;
    set(h1,'Linewidth',2);
    set(h1,'Color',colo);
    h2 = plot(luux,luuy,'k-'); hold on;
    set(h2,'Color',colo);
    h3 = plot(huux,huuy,'k-'); hold on;
    set(h3,'Color',colo);
    %*********
    maxo = max([abs(huux) ; abs(huuy)]) * 1.2;
    axis([-maxo maxo -maxo maxo]);
    plot([-maxo,maxo],[0,0],'k--');
    plot([0,0],[-maxo,maxo],'k--');
    h3 = plot([0,0],[0,maxo],'k-');
    % set(h3,'Linewidth',2);
    %*********
return;

      
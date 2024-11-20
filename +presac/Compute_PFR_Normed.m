function Info = Compute_PFR(Info,Exp,H,plotsteps)
%******* 
%******  function Info = Compute_PFR(Info,Exp,H)
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
       hh = figure('position',[100 100 1800 400]); 
       HA = subplot('position',[0.10 0.15 0.15 0.7]);
       HB = subplot('position',[0.27 0.15 0.15 0.7]);
       HC = subplot('position',[0.45 0.15 0.12 0.7]);
       HD = subplot('position',[0.60 0.15 0.15 0.7]);
       HE = subplot('position',[0.82 0.15 0.15 0.7]);
       %********************
       plot_PFR_unity(Info,HA);   % plot true vs estimate direction
       plot_PFR_angles(Info,HB);  % plot of angle distribution
       plot_PFR_vectors(Info,HC); % plot of PFR vectors
       plot_ProbStat(Info, HD);   % plot of Prob Martix by Targ Freq
       plot_PFR_traces(Info,HE);  % plot the velocity traces
       if (H == 2)  % store result to png and close
          z = getframe(hh);  % current figure
          uname = [Info.pathplot,filesep,'PFR_',Info.tagname,'.png'];
          disp(sprintf('Storing image of PFR graph at %s',uname)); 
          imwrite(z.cdata,uname); % store image for later review
          close(hh);
       end
    else 
       if ~isempty(H)  % H should be a sub-plot handle otherwise
          plot_PFR_unity(Info,H{1});   % plot true vs estimate direction
          plot_PFR_angles(Info,H{2});  % plot of angle distribution
          plot_PFR_vectors(Info,H{3}); % plot of PFR vectors
          plot_ProbStat(Info, H{4});   % plot of Prob Martix by Targ Freq
          plot_PFR_traces(Info,H{5});  % plot the velocity traces
       end
    end
    return;
  end

  % define parameters for doing the PFR analysis
  OpenLoopStart = 0.02;
  OpenLoopEnd = 0.10;
  PostStart = 0.00;
  PostEnd = 0.30;
  PreStart = -0.10;
  PreEnd = 0.00;
  StimSpeed = 15;  % would be better to read from Exp file
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
    openvec = [pxx,pyy];
    
    %******* store vectors, is there a net drift per location
    tg = Info.TargList(k,1); % which of 3 apertures was the target 
    Info.LocationBias{tg} = [Info.LocationBias{tg} ; [pxx,pyy] ];
  end
  for k = 1:NT
      Info.MeanLocBias{k} = median(Info.LocationBias{k});  % use median, reduce outlier
  end
  
  if (plotsteps)
      h = figure(20);
      set(h,'Position',[200 200 (400*NT) 400]);
      for k = 1:NT
          subplot(1,NT,k);   
          for j = 1:size(Info.LocationBias{k},1)
              pxx = Info.LocationBias{k}(j,1);
              pyy = Info.LocationBias{k}(j,2);
              plot([0,pxx],[0,pyy],'k-'); hold on;
          end
          mxx = Info.MeanLocBias{k}(1);
          myy = Info.MeanLocBias{k}(2);
          h = plot([0,mxx],[0,myy],'b-');
          set(h,'Linewidth',2);
          axis tight;
          V = axis;
          maxo = max(abs(V));
          axis([-maxo maxo -maxo maxo]);
          title(sprintf('Location %d',k));
      end
      input('stop');
  end
  
  %***** store velocity traces
  Info.VelTraces = cell(2,2);
  % loop through the list of trials accepted for analyses
  disp('Computing PFR statistics ...');
  Info.PFR_DotList = nan(size(Info.Tlist));  % store dot prods 
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
    
    %******* determine mean drift per target location
    tg = Info.TargList(k,1); % which of 3 apertures was the target
    mxx = Info.MeanLocBias{tg}(1);
    myy = Info.MeanLocBias{tg}(2);  % median PFR shift per location
    vmxx = mxx/(OpenLoopEnd-OpenLoopStart);
    vmyy = myy/(OpenLoopEnd-OpenLoopStart); % median vel vector per loc
    
    %*********
    openvec = [(pxx-mxx),(pyy-myy)];
    
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
    for t = 2:length(plxx)
       velx = (plxx(t) - plxx(t-1))/(pltt(t)-pltt(t-1));  % dva / secs
       vely = (plyy(t) - plyy(t-1))/(pltt(t)-pltt(t-1));  % dva / secs
       vabs = norm([velx,vely]);
       vxx = [vxx velx];
       vyy = [vyy vely];
       vsp = [vsp vabs];
    end
    %***** Apply correction to norm for median drift
    vxx = vxx - vmxx;
    vyy = vyy - vmyy;
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
    %***** Apply correction to norm for median drift
    vxx2 = vxx2 - vmxx;
    vyy2 = vyy2 - vmyy;
    %******* Apply additional smoothing for velocity
    vxx2 = spikestats.gauss_smooth(vxx2,GSIG);
    vyy2 = spikestats.gauss_smooth(vyy2,GSIG);
    vsp2 = spikestats.gauss_smooth(vsp2,GSIG);
   
    %***** go through and NaN away velocities above max speed
    vxx = [vxx2 NaN vxx];
    vyy = [vyy2 NaN vyy];
    vsp = [vsp2 NaN vsp];
    nvxx = vxx;
    nvyy = vyy;
    nvsp = vsp;
    for t = 1:length(nvsp)
        if ~isnan(nvsp(t)) && (nvsp(t) > Max_Eye_Speed)
           ta = max(1,(t-Max_Eye_Gap));
           tb = min(length(nvsp),(t+Max_Eye_Gap));
           %if (ta > 200)
           %    tb = length(nvsp);  %throw out rest of the trial
           % end
           nvxx(ta:tb) = NaN;
           nvyy(ta:tb) = NaN;
        end
    end
    %*******
    
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
    
    %****** LAST is to plot the results on a trial by trial to debug
    if (plotsteps == 1)
        %**********plot some of the data
        h = figure(1); 
        set(h,'position',[100 100 800 800]);  hold off;
        %**** plot the spatial eye trace from fixation
        subplot(2,2,1); hold off;
        plot(xx,yy,'k.-'); hold on;
        h = plot(oxx,oyy,'k.-');
        set(h,'Color',[0.5,0.0,0.0]);  % label open loop trace in gray
        plot(0,0,'r+');  % fixation point
        wid = 8;
        axis([-wid wid -wid wid]);  % bound in -10 to +10 window
        xlabel('X Eye Position');
        ylabel('Y Eye Position');
        title(sprintf('Trial %d, Info %d',tr,k));
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
        
        if (0)
            
            awid = 1.0;
            a1 = startxx - awid;
            a2 = startxx + awid;
            b1 = startyy - awid;
            b2 = startyy + awid;
            plot([a1,a2,a2,a1,a1],[b1,b1,b2,b2,b1],'r-');


            %**** plot the spatial eye trace from fixation
            subplot(2,2,2); hold off;
            plot(xx,yy,'k.-'); hold on;
            %plot the open loop vector
            h = plot(oxx,oyy,'k.-');
            set(h,'Color',[1.0,0.0,0.0]);  % label open loop trace in gray
            % plot(startxx,startyy,'ro');
            % plot([startxx,(startxx+pxx)],[startyy,(startyy+pyy)],'r-');
            plot(0,0,'r+');  % fixation point
            xlabel('X Eye Position');
            ylabel('Y Eye Position');
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
            axis([a1 a2 b1 b2]);  % bound in -10 to +10 window
            h = plot([a1,a2,a2,a1,a1],[b1,b1,b2,b2,b1],'r-');
            set(h,'LineWidth',4);
            title('Zoom in at target');

        else
            %**** plot the time-course of traces 
            subplot(2,2,2); hold off;
            plot(tt,xx,'b.-'); hold on;
            plot(tt,yy,'r.-');
            axis([tt(1) tt(end) -10 10]);
            plot([Tonsac,Tonsac],[-10,10],'k-');
            plot([Toffsac,Toffsac],[-10,10],'k-');
            xlabel('Time (secs)');
            ylabel('Position');

            %plot the target motion and the PFR vector with same origin to compare
            subplot(2,2,3); hold off;
            plot(0,0,'ko'); hold on;
            plot([0,tgvec(1)],[0,tgvec(2)],'b:'); 
            plot([0,openvec(1)],[0,openvec(2)],'r-');
            axis(1.2 * [-endspeed endspeed -endspeed endspeed]);
            xlabel('x velocity');
            ylabel('y velocity');
            title('PFR vector and target motion');

            % aligned the PFR vector along with the motion 
            subplot(2,2,4); hold off;
            plot(0,0,'ko'); hold on;
            plot([0,0],[0,tglength],'b-'); %,'Linewidth',2);
            plot([0,cro],[0,dot],'r-');
            axis(1.2*[-endspeed endspeed -endspeed endspeed]);
            xlabel('x velocity');
            ylabel('y velocity');
            title('PFR vector aligned to motion');
        end
        
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

return;


%***** plot the true motion versus estimated direction on a unity line          
function plot_PFR_unity(Info,H)
  
   subplot(H);
   % pfr directions
   dots = Info.pfrs{1}(:,1); 
   cros = Info.pfrs{1}(:,2); 
   tangs = Info.pfrs{1}(:,4); % actual motion of target apertures
   %*** throw out the outliers (huge PFR due to 2nd sac)
   tglength =  Info.OpenLoopSpeed * (Info.OpenLoopEnd-Info.OpenLoopStart);
   mago = abs( dots + i * cros);
   zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
   %********
   angos = angle( dots(zz) +  i * cros(zz) ) *  (180/pi); %(180/pi);  % go from radians to degs
   tangs = tangs(zz); % match vectors
   %*********
   %*** angos are relative to target, put back into absolute coordinates
   angos = angos + tangs;  % offsets from target angle
   zz = find(angos > 360);  % wrap around circle)
   if ~isempty(zz)
       angos(zz) = angos(zz) - 360;
   end
   zz = find(angos < 0);  % wrap around circle)
   if ~isempty(zz)
       angos(zz) = angos(zz) + 360;
   end
   %*****
   plot(tangs,angos,'b.'); hold on;
   plot([0,360],[0,360],'k-'); 
   axis([0 360 0 360]);
   xlabel('True motion direction');
   ylabel('PFR motion direction');
   title(Info.tagname);
   %****
   
return;

function plot_PFR_angles(Info,H)  % plot of angle distribution
          
   subplot(H);
   %**********
   plot(Info.PFR_TargHistX,Info.PFR_TargHist,'r.-'); hold on;
   plot(Info.PFR_TargHistX,Info.PFR_DistHist,'b.-'); hold on;
   %***********
   xlabel('Angle (degs)');
   ylabel('Counts');
   mederr = Info.PFR_TargErr;
   title(sprintf('Med Targ Err: %5.1f',mederr));
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
  su = nanstd(Info.VelTraces{2,1})/sqrt(size(Info.VelTraces{2,1},1));
  vtt = 1:length(uu);
  h2 = plot(vtt,uu,'r-'); hold on;
  set(h2,'LineWidth',2);
  plot(vtt,(uu+(2*su)),'r-');  % traces plus 2 sem
  plot(vtt,(uu-(2*su)),'r-');  % traces minus 2 sem
  uu = nanmean(Info.VelTraces{2,2});
  su = nanstd(Info.VelTraces{2,2})/sqrt(size(Info.VelTraces{2,2},1));
  vtt = 1:length(uu);
  h2 = plot(vtt,uu,'b-'); hold on;
  set(h2,'LineWidth',2);
  plot(vtt,(uu+(2*su)),'b-');  % traces plus 2 sem
  plot(vtt,(uu-(2*su)),'b-');  % traces minus 2 sem
  plot(vtt,zeros(size(vtt)),'k-');
  axis tight;
  xlabel('Time (ms)');
  ylabel('Velocity Gain');
  
return;

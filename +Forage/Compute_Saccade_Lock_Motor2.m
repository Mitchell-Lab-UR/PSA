function SInfo = Compute_Saccade_Lock_Motor2(Info,Exp,Unit,H,TrialType,Interval)
%******* 
%******  function Info = Compute_DSI(Info,Exp,Unit,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          Unit - spike cluster to run computation on
%***          H - handle to plot window if provided, or if 1 create, else []
%***          TrialType - 1 use BackImage, 2 use Forage (noisetype = 1)
%***
%*** Outputs: Info - it computes the DSI from all trial types,
%***                   and puts confidence intervals on it

  SInfo = [];  % return saccade Info object
  
  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh2 = figure('position',[200 300 1200 300]);
       hh = figure('position',[100 100 400 600]); 
       HA = subplot('position',[0.1 0.35 0.8 0.6]);
       HB = subplot('position',[0.1 0.10 0.8 0.2]);
    else
        if isempty(H)
          return;
        else
          if (length(H) ~= 2)
              disp('Error with subplot pass to Compute_Saccade_Lock_Motor, requires two subplots');
              return;
          else
              HA = H{1};
              HB = H{2};
          end
        end
    end
    %*******
    plot_sacrast(Info,HA);
    plot_sacpsth(Info,HB);
    plot_sacmotor(Info,hh2);
    %********
    if isfloat(H)  && (H == 2)
       z = getframe(hh);  % current figure
       uname = [Info.pathplot,filesep,sprintf('SAC%d_',TrialType),Info.tagname,'.png'];
       disp(sprintf('Storing image of SAC graph at %s',uname)); 
       imwrite(z.cdata,uname); % store image for later review
       close(hh);
    end   
    return;
  end
  
  Info
  Exp
  
  input('Running without plots');
  
  %****** parameters for computing DSI and neuron inclusion
  JackN = 10;             % Jacknife to get confidence intervals
  StimWin = [-0.1,0.2];   % plot stim locked firing rate PSTH 
  Smooth = 5;             % Gaussian sig for smooth PSTH plots
  
  %****** use trial inclusion to build list of trials  *****
  TList = [];
  for fr=1:size(Exp.D,1)
    if (TrialType == 2)  
      if strcmp(Exp.D{fr}.PR.name,'ForageProceduralNoise') 
        if (Exp.D{fr}.PR.noisetype == 1)
           TList = [TList ; fr];             
        end
      end
    else
      if (TrialType == 3)
        %***** any trials where there is almost a gray screen
        if strcmp(Exp.D{fr}.PR.name,'ForageProceduralNoise') 
          if (Exp.D{fr}.PR.noisetype == 3)  % CSD backgrounds
             TList = [TList ; fr];             
          end
        end
        if (strcmp(Exp.D{fr}.PR.name,'InterTrial'))
             TList = [TList ; fr];                 
        end
        %********
        %if strcmp(Exp.D{fr}.PR.name,'BackImage')
        %     if strcmp(Exp.D{fr}.PR.imagefile(13:17),'cloud')  % don't include cloud
        %        TList = [TList ; fr];
        %     end
        % end
        % if strcmp(Exp.D{fr}.PR.name,'ForageStaticLines') 
        %        TList = [TList ; fr];
        % end
        %***********
      else
         if strcmp(Exp.D{fr}.PR.name,'BackImage') 
             if ~strcmp(Exp.D{fr}.PR.imagefile(13:17),'cloud')  % don't include cloud
                TList = [TList ; fr];
             end
         end
      end
    end
  end
  
  %*** take one degree around point (not scaling for RF now)
  %****** get the unit's RF
  if (TrialType == 1) && ~isempty(Info.RFinfo)
     bimage2 = 1920;
     bimage1 = 1080;
     screct = Exp.S.screenRect;  % size of screen (x,y)
     pixperdeg = Exp.S.pixPerDeg;
     bmidx = floor(bimage2/2);
     bmidy = floor(bimage2/2);
     bpixperdegx = (bimage2/screct(3))*pixperdeg;
     bpixperdegy = -(bimage1/screct(4))*pixperdeg;
     bFs = 0.5*(abs(bpixperdegx)+abs(bpixperdegy));
     bsize = floor(0.25*(abs(bpixperdegx)+abs(bpixperdegy)));
     bsize = 2^(ceil(log2(bsize))); % nearest power of 2
     %**** make sure subsequent images match this scale, report if not
     Zbsize = bsize;
     ZbFs = bFs;
     %****
     Gwin = zeros(2*bsize,2*bsize);
     sig2 = (bsize/2)^2;
     for gi = 1:(2*bsize)
         for gj = gi:(2*bsize)
             dist = ((gi-bsize)^2 + (gj-bsize)^2);
             Gwin(gi,gj) = exp(-0.5*(dist/sig2));
             Gwin(gj,gi) = Gwin(gi,gj);
         end
      end
      bw = ((0:bsize)/(2*bsize)) *  bFs;  % cycles per deg 
      %******* can you transform hartley tuning into that coordinate
      %******* system via interpolation?? ****************
      HartU = Info.Hart.uu;
      mu = mean(mean(HartU));
      HartU = (HartU-mu)/mu;
      SpatOris = Info.Hart.SpatOris * (pi/180);
      %*** wrong about Hartley, first is 90 deg, not 0
      if (1)
          SpatOris = SpatOris + (pi/2);
          zz = find(SpatOris >= pi);
          SpatOris(zz) = SpatOris(zz)-pi;
      end
      SpatFrqs = Info.Hart.SpatFrqs;
      bsize2 = bsize/2;
      sigpi2 = (pi/32)^2;
      sigsf2 = (0.5)^2;
      HartW = zeros(bsize,bsize);
      for gi = 1:bsize
         xf = ((gi-(bsize2+1))/(2*bsize)) * bFs;
         for gj = 1:bsize
            yf = ((gj-(bsize2+1))/(2*bsize)) * bFs;
            %****** get distances from RF points
            sf = norm([xf,yf]);
            ango = angle(complex(xf,yf));
            if (ango < 0)
                ango = ango + pi;
            end
            %***************
            sumo = 0;
            wsumo = 0;
            for hi = 1:length(SpatOris)
                dango = ango-SpatOris(hi);
                if (dango > (pi/2))
                    dango = dango - pi;
                end
                if (dango < -(pi/2))
                    dango = dango + pi;
                end
                for hj = 1:length(SpatFrqs)
                    sdis = sf-SpatFrqs(hj);
                    val = exp(-0.5*(sdis^2)/sigsf2) * exp(-0.5*(dango^2)/sigpi2);     
                    sumo = sumo + (val * HartU(hi,hj));
                    wsumo = wsumo + val;
                end
            end
            HartW(gj,gi) = (sumo/wsumo);  %transpose here
         end
      end
      %********************************
  end
  
  %******* then go to collect all the saccade events
  EventList = [];  % list of times in ephys time for time-locking, and 
                   % second column is RT to next saccade
  %*** store saccade vector per event, and let's bin them too
  SacVecList = [];  % [x,y] vector (column 1+2), bin - integer last column
  N_SacDir = 12;  % bins to code direction
  N_SacEcc = 3;  % bins to code eccentricity, power of 2, thus <2, < 4, <8, <16
  for tr = 1:length(TList)
      fr = TList(tr);  % trial number
      estart = Exp.D{fr}.START_EPHYS;
      % MatStart = Exp.D{fr}.eyeData(6,1);  
      %****** READ IN BACKIMAGE IF DOING ROI GRABS PER FIX
      bimage = [];  % empty default
      if (TrialType == 1) && strcmp(Exp.D{fr}.PR.name,'BackImage') && ...
              ~isempty(Info.RFinfo) && ...
               ~isempty(Info.RFinfo.peak) && ...
              ~isempty(Info.Hart)
          %***********
          % disp(Exp.D{fr}.PR.imagefile)
          bimage = imread(Exp.D{fr}.PR.imagefile);
          bimage = mean(bimage,3)/255; % bimage is y by x!, y=0 top
          %***** rescale to desired size *******
          if (size(bimage,1) ~= bimage1) || (size(bimage,2) ~= bimage2)
              bimage = imresize(bimage,[bimage1 bimage2]);
          end
          %****************
          screct = Exp.S.screenRect;  % size of screen (x,y)
          pixperdeg = Exp.S.pixPerDeg;
          bmidx = floor(size(bimage,2)/2);
          bmidy = floor(size(bimage,1)/2);
          bpixperdegx = (size(bimage,2)/screct(3))*pixperdeg;
          bpixperdegy = -(size(bimage,1)/screct(4))*pixperdeg;
          rfx = Info.RFinfo.peak(1);
          rfy = Info.RFinfo.peak(2);
          %*** take one degree around point (not scaling for RF now)
          bFs = 0.5*(abs(bpixperdegx)+abs(bpixperdegy));
          bsize = floor(0.25*(abs(bpixperdegx)+abs(bpixperdegy)));
          bsize = 2^(ceil(log2(bsize))); % nearest power of 2
      end     
      %*******  
      for kt = 1:size(Exp.D{fr}.slist,1)  % list of flagged saccades
         if (kt < size(Exp.D{fr}.slist,1) )
             RT = (Exp.D{fr}.slist(kt+1,2) - Exp.D{fr}.slist(kt,2));
         else
             RT = 1;  % else set as long to next saccade, 1 second
         end
         etime = estart + Exp.D{fr}.slist(kt,1); % sac onset
         otime = estart + Exp.D{fr}.slist(kt,2);  % sac offset
         %******** saccade information
         sta = Exp.D{fr}.slist(kt,4); % integer, sac start
         eta = Exp.D{fr}.slist(kt,5); % integer, sac end
         dx = Exp.D{fr}.eyeSmo(eta,2) - Exp.D{fr}.eyeSmo(sta,2);
         dy = Exp.D{fr}.eyeSmo(eta,3) - Exp.D{fr}.eyeSmo(sta,3);
         mag = norm([dx,dy]);
         ang = (angle(complex(dx,dy))+pi)*(180/pi);
         bin = ( (floor(log2(max(mag,1)))) * N_SacDir);
         bin = bin + 1 + floor(ang/(360/N_SacDir));
         SacVecList = [SacVecList ; [dx,dy,mag,ang,bin]];
         %**************************
         %**************************
         %******* if SacImage, find fixation patch ranking (pref to non-p)
         ImCorr = NaN;
         if (TrialType == 1) && ~isempty(bimage)
             %*** compute the post-saccadic eye position
             endx = Exp.D{fr}.eyeSmo(eta,2);
             endy = Exp.D{fr}.eyeSmo(eta,3);
             rendx = endx + rfx;
             rendy = endy + rfy;
             bx = floor( bmidx + (bpixperdegx * rendx));
             by = floor( bmidy + (bpixperdegy * rendy));
             breakbounds = 0;
             bx1 = (bx-bsize);
             if (bx1 <= 0) 
                 breakbounds = 1;
             end
             bx2 = (bx+bsize-1);
             if (bx2 > size(bimage,2)) 
                 breakbounds = 1;
             end
             by1 = (by-bsize);
             if (by1 <= 0) 
                 breakbounds = 1;
             end
             by2 = (by+bsize-1);
             if (by2 > size(bimage,1)) 
                 breakbounds = 1;
             end
             %***
             %******* show image figure and draw samples
             if ~breakbounds
                 if (0)  % for debugging
                    fp = ((2*pi)/bsize);
                    [X,Y] = meshgrid((1:(2*bsize)),(1:(2*bsize)));
                    imo = 0.5*(1+cos(fp*X));
                    figure(10);
                    imagesc(imo);
                 end
                 imo = bimage(by1:by2,bx1:bx2);
                 imo = 0.5 + ((imo - mean(mean(imo))) .* Gwin);  % window image
                 imu = ( (imo-mean(mean(imo))) )'; % vertical inverted
                 ffimo = abs(fft2(imu));
                 fsize = size(ffimo,1);
                 bsize2 = bsize/2;
                 fimo = zeros(bsize,bsize);
                 fimo(1:bsize2,1:bsize2) = ffimo((fsize-bsize2+1):fsize,(fsize-bsize2+1):fsize);
                 fimo((bsize2+1):bsize,1:bsize2) = ffimo(1:bsize2,(fsize-bsize2+1):fsize);
                 fimo(1:bsize2,(bsize2+1):bsize) = ffimo((fsize-bsize2+1):fsize,1:bsize2);
                 fimo((bsize2+1):bsize,(bsize2+1):bsize) = ffimo(1:bsize2,1:bsize2);
                 % fimo = fimo';
                 %********* note fimo        
                 ImCorr = sum(sum( fimo .* HartW));
          
                 if (0) % debugging
                   hf = figure(10);
                   set(hf,'Position',[100 400 1600 400]);
                   subplot('Position',[0.075 0.15 0.20 0.75]); hold off;
                   colormap('gray');
                   imagesc(bimage,[0 1]); hold on;
                   plot(bmidx,bmidy,'r+');
                   plot([bx1,bx1,bx2,bx2,bx1],[by1,by2,by2,by1,by1],'c-');
                   [rendx,rendy]
                   subplot('Position',[0.300 0.15 0.20 0.75]); hold off;
                   %***** plot the image patch
                   colormap('gray');
                   imagesc(imo,[0,1]); hold on;
                   %*******
                   subplot('Position',[0.525 0.15 0.20 0.75]); hold off;
                   %***** plot the image patch
                   colormap('gray');
                   imagesc(flipud(fimo)); hold on;
                   title(sprintf('ImCorr %f',ImCorr));
                   %*********
                   if (0)
                      subplot('Position',[0.775 0.15 0.20 0.75]); hold off;
                      colormap('gray');
                      imagesc(flipud(fimo));
                   else
                      %****** plot the neurons kernel
                      subplot('Position',[0.775 0.60 0.20 0.35]); hold off;
                      colormap('gray');
                      imagesc(flipud(HartW));
                      %*********
                      subplot('Position',[0.775 0.10 0.20 0.35]); hold off;
                      colormap('gray');
                      imagesc(HartU');
                   end
                   %******
                   bw(1:10)
                   %**** need the spatial freqencies
                   %*******
                   input('check');
                 end
             end
         end
         EventList = [EventList ; [etime,otime,RT,ImCorr]];
         %***************************
      end
      disp(sprintf('Saccade Image matching %d of %d',tr,length(TList)));
  end
  
  
  %****** compute the raster and PSTH ******
  sp = Exp.sp{Unit};
  %*******
  SInfo.EventList = EventList;
  SInfo.StimWin = StimWin;
  SInfo.TrialType = TrialType;
  %****
  %******** resort trials by order of orientation
  StimRast = [];  % keep raster of stim response
  %****** build a column of orienation values (int) and spike rates
  disp('Collecting spike counts per saccade  ...');
  
  spike_times = sp.st;
  event_times = [];
  for k = 1:length(EventList)
      OriInd(k) = EventList(k,2);  % use RT time to sort tirals
      %****** get stim spike count on trial
      onstim = EventList(k,1); 
      staspk = onstim + StimWin(1); % in secs
      finspk = onstim + StimWin(2);
      z = find( (sp.st >= staspk) & (sp.st < finspk) );
      if ~isempty(z)
        difftime = (sp.st(z) - staspk );
        StimRast = [StimRast ; [(difftime+StimWin(1)) (ones(size(difftime))*k)]];
      end
  end 
  
  size(EventList)
  input('sotp');
  
  SInfo.StimRast = StimRast;
  SInfo.OriInd = EventList(:,3);  % trial list of RTs per event 
  % SInfo.OriInd = EventList(:,4);  % trial list by image correlation 
  VM = length(SInfo.OriInd);   % must give PSTH number of trials to fill rast
  [tt,uu,su] = spikestats.CompPSTH(SInfo.StimRast,SInfo.StimWin,VM,Smooth);     
  SInfo.StimTT = tt;
  SInfo.StimUU = uu;   % figure rate before stimulus onset
  SInfo.StimSU = su;  
  %******* compute the base firing rate in pre-stimulus period
  zz = find(tt < 0);  % get baseline firing as time before 0 ms
  if ~isempty(zz)
      SInfo.BaseMu = mean(uu(zz));
  else
      SInfo.BaseMu = NaN;
  end
  %*************
  
  %****************
  if ~isempty(Info.RFinfo)
    %*****  
    EventL = sortrows(EventList,4,'descend'); % sort by RF
    L = length(EventL); 
    %********
    StimRastOff = [];
    StimRastPref = [];
    StimRastNPref = [];
    for k = 1:length(EventL)
      %****** get stim spike count on trial
      onstim = EventL(k,2); % offset
      staspk = onstim + StimWin(1); % in secs
      finspk = onstim + StimWin(2);
      z = find( (sp.st >= staspk) & (sp.st < finspk) );
      if ~isempty(z)
        difftime = (sp.st(z) - staspk );
        StimRastOff = [StimRastOff ; [(difftime+StimWin(1)) (ones(size(difftime))*k)]];
        if (k < floor(0.1*L))
           StimRastNPref = [StimRastNPref ; [(difftime+StimWin(1)) (ones(size(difftime))*k)]];
        end
        if (k > floor(0.9*L))
           kk = (k-floor(0.9*L)); 
           StimRastPref = [StimRastPref ; [(difftime+StimWin(1)) (ones(size(difftime))*kk)]];
        end
      end
    end 
    %******** other computations then
    SInfo.StimRastOff = StimRastOff;
    SInfo.StimRastPref = StimRastPref;
    SInfo.StimRastNPref = StimRastNPref;
  
    %***********
    SInfo.OriInd2 = EventList(:,4);  % trial list by image correlation 
    VM = length(SInfo.OriInd2);   % must give PSTH number of trials to fill rast
    [tt,uu,su] = spikestats.CompPSTH(SInfo.StimRastOff,SInfo.StimWin,VM,Smooth);     
    SInfo.AllStimTT = tt;
    SInfo.AllStimUU = uu;   % figure rate before stimulus onset
    SInfo.AllStimSU = su;  
    %******* top 10 %
    VM = L - floor(0.9*L);
    [att,auu,asu] = spikestats.CompPSTH(SInfo.StimRastPref,SInfo.StimWin,VM,Smooth);     
    SInfo.PrefStimTT = att;
    SInfo.PrefStimUU = auu;   % figure rate before stimulus onset
    SInfo.PrefStimSU = asu;
    %*********
    VM = floor(0.1*L);
    [att,auu,asu] = spikestats.CompPSTH(SInfo.StimRastNPref,SInfo.StimWin,VM,Smooth);     
    SInfo.NPrefStimTT = att;
    SInfo.NPrefStimUU = auu;   % figure rate before stimulus onset
    SInfo.NPrefStimSU = asu;
    %*********
  end
  
  %******** Motor Lock computation ***************
  %%** question, is the peak motor information?  Find the peak, then 
  %*** as for a 50 ms interval centered around the peak
  zz = find( (tt >= -0.05) & (tt < 0.10) );
  tdip = zz( find( uu(zz) == min(uu(zz)) ) );
  dipbin = find( (tt >= (tt(tdip(1))-0.025)) & (tt < (tt(tdip(1))+0.025)) );
  if isempty(Interval)
    %***** constrain, peak must follow dip
    zz = find( (tt >= tt(tdip(1))) & (tt < 0.20) );  % find peak 0 to 200 ms
    if isempty(zz)
      tpeak = tdip(1);
      peakbin = dipbin;
    else
      tpeak = zz( find( uu(zz) == max(uu(zz)) ) );
      peakbin = find( (tt >= (tt(tpeak(1))-0.025)) & (tt < (tt(tpeak(1))+0.025)) );
    end
    SInfo.Interval = [tt(peakbin(1)) tt(peakbin(end))];
  else
    tpeak = mean(Interval);
    peakbin = find( (tt >= Interval(1)) & (tt < Interval(2)) );
  end
  %***
  %******* plot to check it looks right
  if (0)
    figure(10);
    [tt' uu' su']
    plot(tt,uu,'k-'); hold on;
    plot(tt,(uu+(2*su)),'k:');
    plot(tt,(uu-(2*su)),'k:');
    plot([tt(tpeak),tt(tpeak)],[0,max(uu)],'b-');
    plot([tt(tdip),tt(tdip)],[0,max(uu)],'g-');
    input('check');
  end
  %******** now run Sac Events, get spike counts
  SacPeak = cell(N_SacEcc,N_SacDir);
  SacDip = cell(N_SacEcc,N_SacDir);
  for k = 1:length(EventList)
      sbin = SacVecList(k,5);  % bin number
      secc = 1+floor((sbin-1)/N_SacDir);
      if (secc > N_SacEcc)
          continue;  % outside bounds
      else
          sdir = 1+mod((sbin-1),N_SacDir);
          %****** get spike count
          onstim = EventList(k,1); 
          staspk = onstim + tt(peakbin(1)); % in secs
          finspk = onstim + tt(peakbin(end));
          z = find( (sp.st >= staspk) & (sp.st < finspk) );
          pcount = (length(z)/(tt(peakbin(end))-tt(peakbin(1))));
          %*********
          staspk = onstim + tt(dipbin(1)); % in secs
          finspk = onstim + tt(dipbin(end));
          z = find( (sp.st >= staspk) & (sp.st < finspk) );
          dcount = (length(z)/(tt(dipbin(end))-tt(dipbin(1))));
          %********
          SacPeak{secc,sdir} = [SacPeak{secc,sdir} ; pcount];
          SacDip{secc,sdir} = [SacDip{secc,sdir} ; dcount];
          %**************
      end
  end
  
  %********* compute mean and sem
  NSac = zeros(N_SacEcc,N_SacDir);
  UPeak = zeros(N_SacEcc,N_SacDir);
  SPeak = zeros(N_SacEcc,N_SacDir);
  UDip = zeros(N_SacEcc,N_SacDir);
  SDip = zeros(N_SacEcc,N_SacDir);
  %******** values for running anovan
  yy = [];
  yy2 = [];
  g1 = [];
  g2 = [];
  %*****
  for i = 1:N_SacEcc
      for j = 1:N_SacDir
          %******* summarize for peak values
          NSac(i,j) = length(SacPeak{i,j});
          if ~isempty(SacPeak{i,j})
             UPeak(i,j) = nanmean(SacPeak{i,j});
             SPeak(i,j) = nanstd(SacPeak{i,j})/sqrt(length(~isnan(SacPeak{i,j})));
             yy = [yy ; SacPeak{i,j}];
             yy2 = [yy2 ; SacDip{i,j}];
             g1 = [g1 ; (i*ones(size(SacPeak{i,j})))];
             g2 = [g2 ; (j*ones(size(SacPeak{i,j})))];
          else
             UPeak(i,j) = NaN;
             SPeak(i,j) = NaN;
          end
          %****** summarize for dip values
          if ~isempty(SacDip{i,j})
             UDip(i,j) = nanmean(SacDip{i,j});
             SDip(i,j) = nanstd(SacDip{i,j})/sqrt(length(~isnan(SacDip{i,j})));
          else
             UDip(i,j) = NaN;
             SDip(i,j) = NaN;
          end
          %*************
      end
  end
  if (0) % final debugging
    mean(mean(UPeak))
    mean(mean(UDip))
    input('check means');
  end
  %******* PERFORM A TWO-WAY ANOVA ON SAC-MODULATION
  [PeakP] = anovan(yy,{g1,g2},'model','interaction','display','off');
  [DipP] = anovan(yy2,{g1,g2},'model','interaction','display','off');
  SInfo.PeakP = PeakP;
  SInfo.DipP = DipP;
  %********** store results
  SInfo.N_SacEcc = N_SacEcc;
  SInfo.N_SacDir = N_SacDir;
  SInfo.NSac = NSac;
  SInfo.UPeak = UPeak;
  SInfo.SPeak = SPeak;
  SInfo.UDip = UDip;
  SInfo.SDip = SDip;
  SInfo.SacPeak = SacPeak;
  SInfo.SacDip = SacDip;
  SInfo.SacVecList = SacVecList;
  %*********************
return;

%********** Plot motor tuning at peak and dip, z-scored rates
function plot_sacmotor(Info,H)
   figure(H);
   %****
   subplot('position',[0.075 0.15 0.2 0.70]);
   plot(Info.SacVecList(:,1),Info.SacVecList(:,2),'k.','Markersize',2); hold on;
   axis([-16 16 -16 16]);
   %**********
   subplot('position',[0.325 0.15 0.2 0.70]);
   imagesc(Info.NSac); hold on;
   colorbar;
   %*****
   subplot('position',[0.550 0.15 0.2 0.70]);
   upp = nanmean(nanmean(Info.UPeak));
   zpeak = (Info.UPeak - upp) ./ Info.SPeak;
   imagesc(zpeak,[-8 8]); hold on;
   colorbar;
   title(sprintf('Peak(%5.3f,%5.3f,%5.3f)',Info.PeakP));
   %****
   subplot('position',[0.775 0.15 0.20 0.70]);
   udd = nanmean(nanmean(Info.UDip));
   zdip = (Info.UDip - udd) ./ Info.SDip;
   imagesc(zdip,[-8 8]); hold on;
   colorbar;
   title(sprintf('Peak(%5.3f,%5.3f,%5.3f)',Info.DipP));
   %********
return;

%******* below is the plot of stimulus locked PSTH
function plot_sacrast(Info,H)
   %******* This plots 1) firing rate as function of saccade locking,
   %sorted by the RT (RT stored in OriInd)
   subplot(H);
   %****************
   VM = length(Info.OriInd);
   if isempty(Info.StimRast)   % no spikes from unit in the exp
       return;
   end
   %******************
   OriInd2 = flipud(Info.OriInd2);  % flip ascending and descending in raster
   h = spikeplot.PlotTickRaster(Info.StimRastOff,OriInd2,0); hold on;
   set(h,'Linewidth',2);
   %*** superimpose a mean PSTH with error bars into same plot
   tt = Info.AllStimTT;
   uu = Info.AllStimUU;
   su = Info.AllStimSU;
   basemu = Info.BaseMu;
   %***************
   axis([Info.StimWin(1) Info.StimWin(2) 0 VM]);
   plot([0,0],[0,VM],'k-');
   xlabel('Time (secs)');
   ylabel('Trials');
   if (Info.TrialType == 1)
       title('Background Images');
   else
       title('Grating Stimuli');
   end
   %******************
   
return;

function plot_sacpsth(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   %*** superimpose a mean PSTH with error bars into same plot
   tt = Info.AllStimTT;
   uu = Info.AllStimUU;
   su = Info.AllStimSU;
   basemu = Info.BaseMu;
   %********* now plot results
   VM = max( uu + (2*su))*1.1;
   h2 = plot(tt,uu,'k-',tt,(uu+(2*su)),'k-',tt,(uu-(2*su)),'k-'); hold on;
   %*** plot baseline rate line
   zz = find(tt < 0);  % get baseline firing as time before 0 ms
   if ~isnan(basemu)
      vmu = basemu;
      plot([Info.StimWin(1),Info.StimWin(2)],[vmu,vmu],'b--');
      h3 = text((Info.StimWin(1)+0.10*(Info.StimWin(2)-Info.StimWin(1))),...
                 (vmu*1.0),sprintf('%4.1f hz',basemu));
      set(h3,'Fontsize',12);
      set(h3,'Color',[0,0,1]);
      set(h3,'Fontweight','bold');
   end
   %******
   tt = Info.PrefStimTT;
   uu = Info.PrefStimUU;
   su = Info.PrefStimSU;
   h2 = plot(tt,uu,'b-','Linewidth',2);
   h2 = plot(tt,(uu+(2*su)),'b-',tt,(uu-(2*su)),'b-'); hold on;
   VM = max( VM, max( uu + (2*su))*1.1);
   %*********
   tt = Info.NPrefStimTT;
   uu = Info.NPrefStimUU;
   su = Info.NPrefStimSU;
   h2 = plot(tt,uu,'g-','Linewidth',2);
   h2 = plot(tt,(uu+(2*su)),'g-',tt,(uu-(2*su)),'g-'); hold on;
   VM = max( VM, max( uu + (2*su))*1.1);
   %***************
   axis([Info.StimWin(1) Info.StimWin(2) 0 VM]);
   plot([0,0],[0,VM],'k-');
   xlabel('Time (secs)');
   ylabel('Rate (sp/s)');
   %**********
   
return;

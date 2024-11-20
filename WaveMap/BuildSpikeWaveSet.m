function WF = BuildSpikeWaveSet(osp,PHYClu,PlotIt)
%***************
% function WF = reformatTemplateWave(osp,PlotIt)
%*****
%   osp - spike structure output via kilosort and Jake's pipeline
%   PHYClu - if available, label of unit isolations
%   Plotit - if you wish to compare full template unit by unit
%*****
%   WF - reformatted structure without missing channels (replace NaN)
%*****

    %***** find max waveform of each unit 
    [~,max_site] = max(max(abs(osp.temps),[],2),[],3);
    for curr_template = 1:size(osp.temps,1)
	    allwaves(curr_template,:) = ...
	        100*double(osp.temps(curr_template,:,max_site(curr_template)));
    end
    tempPerClu = findTempForEachClu(osp.clu, osp.spikeTemplates);
    tags = 1 + tempPerClu( osp.cids + 1);
    chans = max_site(tags); 
    waves = allwaves(tags,:);  % max wave per cids
    
    %*********** now, for each unit build a CxTxW wave profile
    ChanMin = min(abs(diff(osp.ycoords)));  % inter-electrode distance
    ShankMin = max(abs(diff(osp.xcoords)));  % inter-shank distance
    MaxY = max(osp.ycoords);
    MaxS = max(osp.xcoords);
    NChan = floor(MaxY/ChanMin);
    NShank = floor(MaxS/ShankMin)+1;
    %*******
    disp('  ');
    disp('****************************');
    disp(sprintf('Contact spacing %d  Shank spacing %d',ChanMin,ShankMin));
    disp(sprintf('Channels/ shank %d  Shank number  %d',NChan,NShank));
    disp('****************************');
    disp('Reformatting into standard grid spacing on waveforms');
    %*** build coordinates to chan number without any gaps
    fullchans = [];
    TotChan = (NChan*NShank);
    for k = 1:TotChan
        vk = mod((k-1),NChan);
        hk = floor((k-1)/NChan);
        xo = 0 + (ShankMin * hk);
        yo = MaxY - (vk*ChanMin);
        %****** repaired full space with real channel numbers
        fullchans = [fullchans ; [k xo yo]];
    end
    %******** no go through and rebuild waves with generic structure
    WF.fullchans = fullchans;  % coordinate system of channels
    WF.WaveNum = size(waves,1);
    WF.WaveSamp = size(waves,2);
    WF.NChan = NChan;
    WF.NShank = NShank;
    WF.TotChan = TotChan;
    WF.WaveSet = NaN(WF.WaveNum,WF.TotChan,WF.WaveSamp);
    WF.Iso = NaN(WF.WaveNum,1);
    WF.PeakChan = NaN(WF.WaveNum,1);
    WF.Depth = NaN(WF.WaveNum,1);
    WF.OldDepth = NaN(WF.WaveNum,1);
    WF.Duration = NaN(WF.WaveNum,1);
    WF.EarlyPeak = NaN(WF.WaveNum,1);
    WF.ChanSpan = NaN(WF.WaveNum,1); % number of channels w 50% range as peak
    WF.ISIdt = 0.00025;
    WF.ISItbin = 0:WF.ISIdt:0.100;
    WF.ISIviolations = NaN(WF.WaveNum,1);
    WF.ISItotal = NaN(WF.WaveNum,1);
    WF.ISIdist = NaN(WF.WaveNum,length(WF.ISItbin));
    WF.ISIfit = NaN(WF.WaveNum,5);  % dual fit parameters
    disp('Recomputing normalized waveforms in standard grids');
    %******   
    for uni = 1:WF.WaveNum
       wav = squeeze(double(osp.temps(tags(uni),:,:)));  
       %*** find range of channels 50% of peak
       chansp = range(wav);
       zz = find( chansp > (0.5*max(chansp)) );
       WF.ChanSpan(uni,1) = length(zz);
       %******
       for k = 1:size(wav,2)
          xo = osp.xcoords(k);
          yo = osp.ycoords(k);
          realchan = find( (WF.fullchans(:,2) == xo) & (WF.fullchans(:,3) == yo) );
          WF.WaveSet(uni,realchan,:) = wav(:,k)';
          if (k == chans(uni))
              WF.PeakChan(uni,1) = realchan;  % channel of peak waveform
              WF.OldDepth(uni,1) = osp.clusterDepths(1,uni); 
              WF.Depth(uni,1) = WF.fullchans(realchan,3);  % y coordinate of peak channel
              waveform = wav(:,k)';
              ztrough = find( waveform == min(waveform) );
              if ~isempty(ztrough) 
                  bzz = ztrough(1):length(waveform);
                  zlatepeak = find( waveform(bzz) == max(waveform(bzz)) );
                  WF.Duration(uni,1) = zlatepeak(1)*(1000/30); % in micro secs, peak after trough
                  latepeak = waveform(bzz(zlatepeak));
                  bzz = 1:(ztrough(1)-1);
                  zearlypeak = find( waveform(bzz) == max(waveform(bzz)) );
                  earlypeak = waveform(bzz(zearlypeak));
                  WF.EarlyPeak(uni,1) = (earlypeak - latepeak)/(abs(earlypeak)+abs(latepeak));
              end
          end
       end 
       %****************
       if ~isempty(PHYClu)
          if (size(PHYClu,1) == (WF.WaveNum+1))  % no extra unit in table
              uu = uni+1;
          else
              if (size(PHYClu,1) == WF.WaveNum)  % no extra unit in table
                  uu = uni;
              else
                  if (size(PHYClu,1) == (WF.WaveNum-1))
                      uu = uni;
                  else
                      uu = uni; 
                  end
              end
          end
          if ( uu > size(PHYClu,1))
             disp(sprintf('UNABLE TO RECOVER ISOLATION, %d',uni));
             sp = [NaN NaN];
          else
             if strcmp(PHYClu.group(uu),'good')  %% what to set this as? 
               WF.Iso(uni,1) = 4; 
             end
             if strcmp(PHYClu.group(uu),'mua')  %% what to set this as? 
               WF.Iso(uni,1) = 1;
             end
             if strcmp(PHYClu.group(uu),'noise')  %% what to set this as? 
               WF.Iso(uni,1) = 0;
             end
             %***** found label ok, get spikes for that unit
             clusto = int32( PHYClu.cluster_id(uu) );
             zsp = find( osp.clu == clusto );
             sp = osp.st(zsp);  % spike times (in secs)
          end
       end
       %*************
       % compute the ISI distrib for the cluster
       isi = hist(diff(sp),WF.ISItbin);
       WF.ISItotal(uni) = sum(isi);
       WF.ISIdist(uni,:) = isi;
       %******* for ISI violations, look at violations within 25 ms
       zz = find( WF.ISItbin <= 0.0010);
       zzt = find( WF.ISItbin < WF.ISItbin(end) );  % don't including ends
       WF.ISIviolation(uni) = 100*(sum(isi(zz))/sum(isi(zzt)));
       %***** Fit ISI process model to get refractory period and burstiness
       if (WF.Iso(uni,1) > 1)
         TP = floor(1/(1000*WF.ISIdt));
         TP
         isidist = squeeze(WF.ISIdist(uni,:));
         isidist = isidist(1:(end-1));  % throw out end (has many points)
         isidist = isidist / sum(isidist);  % norm to frequenzy
         finfo = Fit_Wave_Isi(isidist,TP,1);
         finfo
         finfo.zall
         finfo.all
         input('check');
       end
       %**********
    end
    disp('WaveSet computed, done!');
    %********** Resort the channels by order before storage! ***
    if (0)
      disp('Re-ordering so waves go by channel order');
      chanvals = [(1:WF.WaveNum)' WF.PeakChan];
      chansort = sortrows(chanvals,2);      
      RR = chansort(:,1);  % redo the order
      WF.WaveSet(:,:,:) = WF.WaveSet(RR,:,:);
      WF.Iso = WF.Iso(RR);
      WF.PeakChan = WF.PeakChan(RR);
      chans = chans(RR);
      waves = waves(RR,:,:);
    end
    %*************
    disp('****************************');
    
    %************
    if (PlotIt)
        StepPlotWaveSet(WF,osp);
    end
    %***********
    % Notes on format:  1 starts at electrode tip (deep) to 32 (shallow)
    %                  33 starts tip of second (deep) to 64 (shallow)
    %**** this will be important when matching to CSD analyses
    %***********************
    
return;
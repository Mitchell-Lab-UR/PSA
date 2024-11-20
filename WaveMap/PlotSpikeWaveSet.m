function PlotSpikeWaveSet(WF)
%***************
%** function PlotSpikeWaveSet(WF)
%**************
%   WF - struct with information about waveforms across electrode
%   WF.fullchans - three columns: chan number, x coord, y coord
%   WF.WaveNum   - total number of waveforms
%   WF.WaveSamp  - size of waveform (samples at 30 khz)
%   WF.NChan - number of channels per shank
%   WF.NShank - number of shanks
%   WF.TotChan - total channels (dead channels filled in, if any)
%   WF.WaveSet - size (WF.WaveNum,WF.TotChan,WF.WaveSamp) - waveforms
%   WF.Iso - isolation of waveform (WF.WaveNum,1) with (4-good, 1-mua, 0-hash)
%   WF.PeakChan - channel with peak amp waveform (WF.WaveNum,1)
%   WF.Depth - depth of the peak channel (from peak only)
%   WF.OldDepth - average depth from multiple channel peaks
%   WF.Duration - depth of waveform at channel peak (peak dur after trough)
%   WF.EarlyPeak - positive value (-1 to 1) evidence of peak before trough
%   WF.ChanSpan - number of channels with range within 50% of peak
%   WF.ISIdist - (WF.Wavenum,N_ISI) ISI distribution within 25 ms
%   WF.ISIdt = 0.00025;  % delta t in getting ISI dist
%   WF.ISItbin = time bins in ISI distributions
%   WF.ISIviolations = (WF.Wavenum,1) percent of ISI under 1 ms relative to 25 ms
%   WF.ISItotal = (WF.Wavenum,1) total number of spikes  
%%

  chanum = WF.TotChan;
  shanknum = WF.NChan;
  wavesize= WF.WaveSamp;
  %********
  if (1)  % filter to a subset of units
     zz = find( WF.Iso > 1);  % take on isolated units
     waves = WF.WaveSet(zz,:,:);
     channel = WF.PeakChan(zz);
     unitnum = length(zz); 
  else
     waves = WF.WaveSet;
     channel = WF.PeakChan;
     unitnum = WF.WaveNum;
  end
  
  %*** get order of channels
  chanvals = [(1:unitnum)' channel];
  chansort = sortrows(chanvals,2);
  
  %Plot each unit 
  hf = figure;
  set(hf,'Position',[100 100 1500 750]);
  for i=1:unitnum
    ii = chansort(i,1);
    ipeak = chansort(i,2);
    colo = rand(1,3);
    if (0) % plot peak channel
      chan = mod((ipeak-1), shanknum);  % 0 to (shanknum-1)
      wavePos = ((chan*unitnum)+i);
      mywaves = squeeze(waves(ii,ipeak,:));    
      subplot(shanknum,unitnum,wavePos);
      plot(mywaves, 'Color',colo);
      set(gca,'Visible','off')
      set(gcf,'color','w');
    else  % plot over 32 channels
      mywaves = squeeze(waves(ii,ipeak,:));    
      maxo = max(mywaves);
      mino = min(mywaves);
      if (0)
         vset = 1:shanknum;
         thresh = 0.1;
      else
         ipeakmod = mod((ipeak-1),shanknum)+1;
         vlow = max((ipeakmod-3),1);
         vhigh = min((ipeakmod+3),shanknum);
         vset = vlow:vhigh;
         thresh = 0.0;
      end
      %**********
      for vi = vset
         chan = vi-1;  % 0 to (shanknum-1)
         wavePos = ((chan*unitnum)+i);
         if (ipeak > 32)
             ip = 32+vi;
         else
             ip = vi;
         end
         mywaves = squeeze(waves(ii,ip,:));
         imaxo = max(mywaves);
         imino = min(mywaves);
         if ((imaxo-imino) > (thresh*(maxo-mino)))
           subplot(shanknum,unitnum,wavePos);
           plot(mywaves, 'Color',colo); hold on;
           if (vi == 1)
               plot([0,wavesize],[maxo,maxo],'k-','LineWidth',2);
           end
           if (vi == shanknum)
               plot([0,wavesize],[mino,mino],'k-','LineWidth',2);
           end
           axis([0 wavesize mino maxo]);
           set(gca,'Visible','off')
           set(gcf,'color','w');
         end
      end
      if (1)  % replot the peak
        chan = mod((ipeak-1), shanknum);  % 0 to (shanknum-1)
        wavePos = ((chan*unitnum)+i);
        mywaves = squeeze(waves(ii,ipeak,:));    
        subplot(shanknum,unitnum,wavePos);
        plot(mywaves,'Color',colo,'LineWidth',1.5);
        axis([0 wavesize mino maxo]);
        set(gca,'Visible','off')
        set(gcf,'color','w');
      end
      %******
      disp(sprintf('Plotting unit %d',i));
      % input('check');
    end
  end
  
end


function StepPlotWaveSet(WF,osp)
%*** function StepPlotWaveSet(WF)
%**
%**  WF - WaveformSet struct (info of spike waveforms derived
%**       from osp struct and Exp.PHYClu structs on isolation quality
%**  osp - original osp struct, for comparison of template data
%**
%** Goes through unit by unit to plot waveform across channels
%** and report stats (ISI dist) on isolation quality
%**
        %****** need to recompute tags
        if ~isempty(osp)
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
        end

        hf = figure;
        set(hf,'Position',[100 100 600 900]);
        for uni = 1:WF.WaveNum
           %***** plot original description of waves
           if ~isempty(osp)
               subplot('Position',[0.1 0.05 0.35 0.70]); 
               wav = 100*squeeze(double(osp.temps(tags(uni),:,:)));
               tt = 1:size(wav,1);
               hold off;
               for k = 1:size(wav,2)
                xo = osp.xcoords(k);
                yo = osp.ycoords(k);
                plot((xo+tt),(yo+wav(:,k)'),'k-'); hold on;
               end
               xo = osp.xcoords(chans(uni));
               yo = osp.ycoords(chans(uni));
               plot((xo+tt),(yo+waves(uni,:)),'r-'); % red on peak chan
               text(100,1000,sprintf('Ch %d',chans(uni)),'Fontsize',12);
           end
           %**********
           %******* OK, now plot waves from the normalized format with no deads
           subplot('Position',[0.60 0.05 0.35 0.70]);
           wav = 100*squeeze(double(WF.WaveSet(uni,:,:)));
           tt = 1:WF.WaveSamp;
           hold off;
           for k = 1:WF.TotChan
              xo = WF.fullchans(k,2);
              yo = WF.fullchans(k,3);
              plot((xo+tt),(yo+100*squeeze(WF.WaveSet(uni,k,:))),'k-'); hold on;
           end
           xo = WF.fullchans(WF.PeakChan(uni),2);
           yo = WF.fullchans(WF.PeakChan(uni),3);
           plot((xo+tt),(yo+100*squeeze(WF.WaveSet(uni,WF.PeakChan(uni),:))),'r-'); % red on peak chan
           text(100,1000,sprintf('Ch %d',WF.PeakChan(uni)),'Fontsize',12);
           colo = [0,0.7,0.7];
           text(110,900,sprintf('Iso %d',WF.Iso(uni,1)),'Color',colo);
           text(90,850,sprintf('Dep %5.1f',WF.Depth(uni,1)),'Color',colo);
           text(90,800,sprintf('ODp %5.1f',WF.OldDepth(uni,1)),'Color',colo);
           text(100,750,sprintf('Dur %5.1f',WF.Duration(uni,1)),'Color',colo);
           text(100,700,sprintf('ErPk %4.2f',WF.EarlyPeak(uni,1)),'Color',colo);
           text(100,650,sprintf('ChSp %d',WF.ChanSpan(uni,1)),'Color',colo);
           %**** plot the ISI distribution for this unit (ISO looks right?)
           subplot('Position',[0.10 0.825 0.35 0.125]); hold off;
           % show the ISI distrib for the cluster
           vx = WF.ISItbin;
           isi = WF.ISIdist(uni,:);
           vx = vx(1:(end-1));
           isi = isi(1:(end-1));
           isi = (100 * (isi / sum(isi)));
           bar(1000*vx,isi,1); hold on;
           V = axis;
           plot([1,1],[0,V(4)],'r-');
           xlabel('Time (ms)');
           ylabel('Counts');
           title(sprintf('Viol(%4.2f)  N(%d)',WF.ISIviolation(uni),WF.ISItotal(uni)));
           %******** find the strongest peak cross-correlation?
           subplot('Position',[0.60 0.825 0.35 0.125]); hold off;
           xlabel('Time (ms)');
           ylabel('Cross-corr');   % could plot within 1ms cross-corr per units
                                   % norming by the shuffle-corrector (1
                                   % so a line a 1 is default, if high peak
                                   % is shown away from that line, violate)
           %*********
           input('check');
        end


return;
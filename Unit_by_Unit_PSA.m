MTC=0;

if (1) % Jude's Laptop
    % I think this should be right for Jude's laptop
    StorePath = 'C:\Users\amybu\Box\MTC\Info_MT';  % MT data
        StorePath2 = 'C:\Users\amybu\Box\MTC\Info_MTC';  % MTC data

else
        StorePath = 'C:\Users\amybu\Box\V1Hart\Info_MT2';  % MT data
        StorePath2 = 'C:\Users\amybu\Box\V1Hart\Info_MTC';  % MTC data
        %StorePath = 'C:\Users\amybu\Box\V1Hart\Info_Foveal';  % Foveal data
end


InfoNames = dir([StorePath,filesep,'*']); %MT
InfoNames2 = dir([StorePath2,filesep,'*']);

%Run Unit by Unit
MTC=1;
if MTC==1
    names=InfoNames2(3:end);
else
    names=InfoNames(3:end);
end

%%
for k = 1:length(names)
          a_rand = names(randperm(length(names)));
          fname = a_rand.name;
          disp(sprintf('Loading unit %s',fname));
          load([StorePath2,filesep,fname]);   % load Info struct

            

%end


%% 
    if (Info.Included)
      hf = figure(20); 
      set(hf,'Position',[100 100 1200 800]);
      
      %vmax = 160;
      %vmax2 = 160;
      
      %****** Plot parameters of the timing
      if (1)
          sta = -0.05;
          stn = 0.25; % 0.30;
          bd = 0.20; %0.25;   
          bd2 = 0.0;
          asta = -0.025;
          astn = 0.125; %0.125;
          bsta = -0.05; %-0.075;
          bstn = 0.50;
      else
          sta = -0.05;
          stn = 0.30;
          bd = 0.25;   
          bd2 = 0.0;
          asta = -0.025;
          astn = 0.125;
          bsta = -0.075;
          bstn = 0.50;          
      end
      convo = 0;
      if (0) % show deconvolution
          convo = 1;
          astn = 0.200;  
          bsta = -0.250; 
          stn = 0.50;
          bd = 0.45;
      end

%%      %***** Plot Raster with Saccade Locked timing
      if (1) %exist('Info.StimOn', 'var')
          
         subplot('Position',[0.10 0.40 0.35 0.50]);
         rasto = cell(2,2);
         oro = cell(2,2);
         for ik = 1:2
            dn = 1:floor(length(Info.StimOn.OriIndRT{ik})/Info.MatchN);
            for k = 1:length(dn)
             rr = find( (Info.StimOn.Rast{ik}(:,2) == dn(k)) & ...
                        (Info.StimOn.Rast{ik}(:,1) >= asta) & ...
                        (Info.StimOn.Rast{ik}(:,1) < astn) );
             rasto{1,ik} = [rasto{1,ik} ; Info.StimOn.Rast{ik}(rr,:)];
             oro{1,ik} = [oro{1,ik} ; Info.StimOn.OriIndRT{ik}(dn(k))];
            end
            rasto{1,ik}(:,1) = rasto{1,ik}(:,1) + bd2;
         end
         for ik = 1:2
            dn = 1:floor(length(Info.SacOn.OriIndRT{ik})/Info.MatchN);
            for k = 1:length(dn)
             rr = find( (Info.SacOn.Rast{ik}(:,2) == dn(k)) & ...
                        (Info.SacOn.Rast{ik}(:,1) >= bsta) & ...
                        (Info.SacOn.Rast{ik}(:,1) < bstn) );
             rasto{2,ik} = [rasto{2,ik} ; Info.SacOn.Rast{ik}(rr,:)];
             oro{2,ik} = [oro{2,ik} ; Info.SacOn.OriIndRT{ik}(dn(k))];
            end
            rasto{2,ik}(:,1) = rasto{2,ik}(:,1) + bd;
         end
         maxtrials = length(oro{1,1}) + length(oro{1,2});
         %*****
         aa = [sta stn stn sta sta];
         bb = [maxtrials maxtrials 0 0 maxtrials];
         fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
         %*** add spike rasters here with other sac timings
         h2 = spikeplot.PlotTickRaster(rasto{1,2},oro{1,2},0); hold on;
         VM = length(Info.SacOn.OriInd{2});
         MN = Info.MatchN;
         [tt12,uu12,su12] = spikestats.CompPSTH_Boot(Info.SacOn.Rast{2},Info.SacOn.LockWin,VM,MN,Info.SacOn.Smooth);     
         set(h2,'Linewidth',1);
         set(h2,'Color',[0,0,1]);
         voffset = length(oro{1,2});
         %*********
         h3 = spikeplot.PlotTickRaster(rasto{1,1},oro{1,1},voffset);
         VM = length(Info.SacOn.OriInd{1});
         [tt11,uu11,su11] = spikestats.CompPSTH_Boot(Info.SacOn.Rast{1},Info.SacOn.LockWin,VM,MN,Info.SacOn.Smooth);     
         set(h3,'Linewidth',1);
         set(h3,'Color',[1,0,0]);
         %   plot(oro{1,2},(voffset+1:length(oro{1,2})),'k.','Markersize',10);
         %********
         h4 = spikeplot.PlotTickRaster(rasto{2,2},oro{2,2},0); hold on;
         VM = length(Info.StimOn.OriInd{2});
         [tt22,uu22,su22] = spikestats.CompPSTH_Boot(Info.StimOn.Rast{2},Info.StimOn.LockWin,VM,MN,Info.StimOn.Smooth); 
         set(h4,'Linewidth',1);
         set(h4,'Color',[0,0,1]);
         % plot(oro{1,1},1:length(oro{1,1}),'k.','Markersize',10);
         voffset = length(oro{2,2});
         %*********
         h5 = spikeplot.PlotTickRaster(rasto{2,1},oro{2,1},voffset);
         VM = length(Info.StimOn.OriInd{1});
         [tt21,uu21,su21] = spikestats.CompPSTH_Boot(Info.StimOn.Rast{1},Info.StimOn.LockWin,VM,MN,Info.StimOn.Smooth);     
         set(h5,'Linewidth',1);
         set(h5,'Color',[1,0,0]);
  
         %***************
         axis([sta stn 0 maxtrials]);
         axis off;
         %*******
         aa = [0.04 0.09 0.09 0.04 0.04];
         bb = [maxtrials maxtrials 0 0 maxtrials];
         fill(aa,bb,[0.5,0.5,0.5],'Linestyle','none','FaceAlpha',0.15); hold on;
         aa = [(bd-0.03) (bd+0.03) (bd+0.03) (bd-0.03) (bd-0.03)];
         bb = [maxtrials maxtrials 0 0 maxtrials];
         fill(aa,bb,[0.5,0.5,0.5],'Linestyle','none','FaceAlpha',0.15); hold on;
         plot([0,0],[0,maxtrials],'k-');
         text(-0.04,maxtrials*1.10,'Stimulus','Fontsize',18);
         text(-0.03,maxtrials*1.05,'Onset','Fontsize',18);
         plot([bd,bd],[0,maxtrials],'k-');
         text(bd-0.03,maxtrials*1.10,'Saccade','Fontsize',18);
         text(bd-0.02,maxtrials*1.05,'Onset','Fontsize',18);
         plot([-0.04,-0.04],[0,100],'k-','Linewidth',2);
         text(-0.06,5,'100 trials','Rotation',90,'Fontsize',18);
         %********
      end
      
      %****** plot the PSTH of single unit below
      if (1)
         subplot('Position',[0.10 0.075 0.35 0.30]);
         vmax = Info.BaseMu *10;
         %******
%          tt1 = Info.SacOn.LockTT{1};
%          tt2 = Info.StimOn.LockTT{1};
%          att = cell(2,2);
%          att{1,1} = Info.SacOn.LockUU{1};
%          att{1,2} = Info.SacOn.LockSU{1};
%          att{2,1} = Info.StimOn.LockUU{1};
%          att{2,2} = Info.StimOn.LockSU{1};
%          utt = cell(2,2);
%          utt{1,1} = Info.SacOn.LockUU{2};
%          utt{1,2} = Info.SacOn.LockSU{2};
%          utt{2,1} = Info.StimOn.LockUU{2};
%          utt{2,2} = Info.StimOn.LockSU{2};
         %*****
         tt1 = Info.SacOn.LockTT{1};
         tt2 = Info.StimOn.LockTT{1};
         att = cell(2,2);
         att{1,1} = uu11; % Info.SacOn.LockUU{1};
         att{1,2} = su11;
         att{2,1} = uu21; % Info.StimOn.LockUU{1};
         att{2,2} = su21;
         utt = cell(2,2);
         utt{1,1} = uu12; % Info.SacOn.LockUU{2};
         utt{1,2} = su12; % Info.SacOn.LockSU{2};
         utt{2,1} = uu22; % Info.StimOn.LockUU{2};
         utt{2,2} = su22;
         %*****
     
         %*********
         SInt = Info.SacOn.LockInt;
         zz = find( (tt1 >= SInt(1)) & (tt1 < SInt(2)) );
         aa = mean(uu11(zz));
         sa = mean(su11(zz))/sqrt((SInt(2)-SInt(1))*1000/(4*Info.SacOn.Smooth));
         uu = mean(uu12(zz));
         su = mean(su12(zz))/sqrt((SInt(2)-SInt(1))*1000/(4*Info.SacOn.Smooth));
         tval = (aa-uu)/sqrt( sa.^2 + su.^2 );
         pval = 2*(1-normcdf(abs(tval)));
         disp('  ');
         disp(sprintf('PSTH Rate Windows:'));
         disp(sprintf('Sac Onset %5.3f to %5.3f',SInt(1),SInt(2)));
         disp(sprintf('Rate (%4.2f,%4.2f) p(%10.8f)',aa,uu,pval));
         disp(' '); 
         SInt = Info.StimOn.LockInt;
         zz = find( (tt2 >= SInt(1)) & (tt2 < SInt(2)) );
         aa = mean(uu21(zz));
         sa = mean(su21(zz))/sqrt((SInt(2)-SInt(1))*1000/(4*Info.StimOn.Smooth));
         uu = mean(uu22(zz));
         su = mean(su22(zz))/sqrt((SInt(2)-SInt(1))*1000/(4*Info.StimOn.Smooth));
         tval = (aa-uu)/sqrt( sa.^2 + su.^2 );
         pval = 2*(1-normcdf(abs(tval)));
         disp(sprintf('PSTH Rate Windows:'));
         disp(sprintf('Stim Onset %5.3f to %5.3f',SInt(1),SInt(2)));
         disp(sprintf('Rate (%4.2f,%4.2f) p(%10.8f)',aa,uu,pval));
         disp(' '); 
         %********
         
         aa = [sta stn stn sta sta];
         bb = [vmax vmax 0 0 vmax];
         fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
         kk = 1;
         zz = find( (tt1 >= bsta) & (tt1 < bstn) );
         aa = [(bd+tt1(zz)) fliplr(bd+tt1(zz))];
         bb = [(att{kk,1}(zz)+2*att{kk,2}(zz)) fliplr(att{kk,1}(zz)-2*att{kk,2}(zz))];
         fill(aa,bb,[1,0,0],'Linestyle','none','FaceAlpha',0.3); hold on;
         plot(bd+tt1(zz),att{kk,1}(zz),'r-','Linewidth',2); hold on;
         bb = [(utt{kk,1}(zz)+2*utt{kk,2}(zz)) fliplr(utt{kk,1}(zz)-2*utt{kk,2}(zz))];
         fill(aa,bb,[0,0,1],'Linestyle','none','FaceAlpha',0.3); hold on;
         plot(bd+tt1(zz),utt{kk,1}(zz),'b-','Linewidth',2); hold on;
         %**** show time window of analyses
         aa = [0.04 0.09 0.09 0.04 0.04];
         bb = [vmax vmax 0 0 vmax];
         fill(aa,bb,[0.5,0.5,0.5],'Linestyle','none','FaceAlpha',0.15); hold on;
         %******
         kk = 2;
         zz = find( (tt2 >= asta) & (tt2 < astn) );
         aa = [(bd2+tt2(zz)) fliplr(bd2+tt2(zz))];
         bb = [(att{kk,1}(zz)+2*att{kk,2}(zz)) fliplr(att{kk,1}(zz)-2*att{kk,2}(zz))];
         fill(aa,bb,[1,0,0],'Linestyle','none','FaceAlpha',0.3); hold on;
         plot(bd2+tt2(zz),att{kk,1}(zz),'r-','Linewidth',2); hold on;
         bb = [(utt{kk,1}(zz)+2*utt{kk,2}(zz)) fliplr(utt{kk,1}(zz)-2*utt{kk,2}(zz))];
         fill(aa,bb,[0,0,1],'Linestyle','none','FaceAlpha',0.3); hold on;
         plot(bd2+tt2(zz),utt{kk,1}(zz),'b-','Linewidth',2); hold on;
         axis tight;
         %V = axis;
         %axis([-0.05 V(2) 0 vmax]);
         axis off;
         plot([0.002,0.102],[0,0],'k-','Linewidth',2);
         text(0.023,(-0.075*vmax),'100 ms','Fontsize',18);
         plot([0,0],[0.01*vmax,1.0*vmax],'k-');
         % text(-0.04,vmax*0.95,'Stimulus','Fontsize',18);
         % text(-0.03,vmax*0.85,'Onset','Fontsize',18);
         plot([bd,bd],[0.01*vmax,1.0*vmax],'k-');
         % text(bd-0.03,vmax*0.95,'Saccade','Fontsize',18);
         % text(bd-0.02,vmax*0.85,'Onset','Fontsize',18);
         %plot([-0.04,-0.04],[0,50],'k-','Linewidth',2);
         %text(-0.06,5,'50 sp/s','Rotation',90,'Fontsize',18);
         text(-0.09,vmax*0.2,'Mean Rate','Rotation',90,'Fontsize',18);
         %**** show time window of analyses
         aa = [(bd-0.03) (bd+0.03) (bd+0.03) (bd-0.03) (bd-0.03)];
         bb = [vmax vmax 0 0 vmax];
         fill(aa,bb,[0.5,0.5,0.5],'Linestyle','none','FaceAlpha',0.15); hold on;
         %********
      end
      
      %************************
      NPSMOOTH = 1;
      centerwid = 0;
      if (1) %Motion Tuning 
         afit1 = Info.SacOn.OriTune{1};
         ovals = Info.SacOn.OriTune{1}.OriVals;
         ovals = (ovals+22.5)/22.5;
         afit2 = Info.SacOn.OriTune{2};
         LO = length(ovals);
         %****** find prefs
         % smooth on the circle ******
         aa = circsmooth(afit1.Mu,NPSMOOTH);
         uu = circsmooth(afit2.Mu,NPSMOOTH);
         sa = circsmooth(afit1.Sem,NPSMOOTH);
         su = circsmooth(afit2.Sem,NPSMOOTH);
         tot = aa + uu;
         theta = 0:(LO-1);
         theta = theta*(2*pi)/LO;
         wvec = sum( tot .* exp(j*theta) );
         wvec = wvec/abs(wvec);
         %********
         if (centerwid) 
              awvec = sum( aa .* exp(j*theta) );
              awvec = awvec/abs(awvec);
              uwvec = sum( uu .* exp(j*theta) );
              uwvec = uwvec/abs(uwvec);
              %*******
              wvec = 0.5*(awvec+uwvec);
              wvec = wvec / abs(wvec);
         end
         %********
         Wvec2 = wvec;
         ango = angle(wvec);
         if (ango < 0)
             ango = ango + 2*pi;
         end
         wvc = (LO*ango/(2*pi)) + 1;
         wshift = 14;
         wvc = wvc + wshift;
         rval = [];
         for ii = 1:length(ovals)
                    ang = ((ii-1)*2*pi/length(ovals));
                    dot = (real(Wvec2)*cos(ang) + imag(Wvec2)*sin(ang));
                    rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
         end
         yrvals = sortrows(rval,2);
         rord = yrvals(:,1);
         NPSMOOTH = 1;
         ratt = aa(rord);  %smoothnp(afit1.Mu(rord),NPSMOOTH);
         sratt = sa(rord); %smoothnp(afit1.Sem(rord),NPSMOOTH);
         rutt = uu(rord);  %smoothnp(afit2.Mu(rord),NPSMOOTH);
         srutt = su(rord); %smoothnp(afit2.Sem(rord),NPSMOOTH);
         [ampw,umpw,abasew,ubasew,aw,uw] = compute_NP_params(ratt,rutt);
         %***********       
         vmax = 1.1*max([max(afit1.Mu),max(afit2.Mu)]);
         aa = [0 (LO+1) (LO+1) 0 0];
         bb = [vmax vmax 0 0 vmax];

         %******get Von Mises fit params and variances
         fit1 = Info.SacOn.OriTune{1}.Afit;  % constrained same pref
         fit2 = Info.SacOn.OriTune{2}.Afit;  % constrained same pref
         disp('Von Mises Fit Differences');
         diffs = fit1.mu - fit2.mu;
         sems = sqrt( (fit1.sem .^2) + (fit2.sem .^2));
         tvals = diffs ./ sems;
         pvals = normcdf(abs(tvals));
         pvals = 2*(1-pvals);
         disp(sprintf('Von Mises Stats:'));
         disp(sprintf('Base(%4.2f,%4.2f) p(%10.8f)',fit1.mu(1),fit2.mu(1),pvals(1)));
         disp(sprintf('Amp (%4.2f,%4.2f) p(%10.8f)',fit1.mu(2),fit2.mu(2),pvals(2)));
         disp(sprintf('Kappa (%4.2f,%4.2f) p(%10.8f)',fit1.mu(3),fit2.mu(3),pvals(3)));
         disp(' '); 

         %******* plot original tuning
             subplot('Position',[0.55 0.625 0.15 0.30]);

         %******

         fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
         axis([0 LO+1 0 vmax]);
         axis off;
         xx = 1:LO;
         if (wshift > 0)
            xx = [(LO-wshift+1):LO, 1:(LO-wshift)]; % recenter
         end
         %******
         aa = circsmooth(afit1.Mu,NPSMOOTH);
         uu = circsmooth(afit2.Mu,NPSMOOTH);
         sa = circsmooth(afit1.Sem,NPSMOOTH);
         su = circsmooth(afit2.Sem,NPSMOOTH);
         %******
         zaa = [1:LO fliplr(1:LO)];
         zbb = [(aa(xx)+2*sa(xx)) fliplr(aa(xx)-2*sa(xx))];
         fill(zaa,zbb,[1,0,0],'FaceAlpha',0.1,'Linestyle','none'); hold on;
         zaa = [1:LO fliplr(1:LO)];
         zbb = [(uu(xx)+2*su(xx)) fliplr(uu(xx)-2*su(xx))];
         fill(zaa,zbb,[0,0,1],'FaceAlpha',0.1,'Linestyle','none'); hold on;
         plot(ovals,aa(xx),'r.-','Markersize',10); hold on;
         plot(ovals,uu(xx),'b.-','Markersize',10); hold on;
         plot([wvc,wvc],[0,vmax],'k--','Linewidth',1);
         text(4,-0.1*vmax,'Motion direction','Fontsize',14);
         plot([0.0,0.0],[0,30],'k-','Linewidth',2);
         text(-1.0,3,'30','Fontsize',12,'Rotation',90);
         text(-2.5,0.1*vmax,'Rate (sp/s)','Rotation',90,'Fontsize',14);
         %*****
         text(4,1.15*vmax,'Motion Tuning','Fontsize',16);
         text(3,1.05*vmax,'at Saccade Onset','Fontsize',16);
         %text(-0,1.15*vmax,'Non-Parametric(NP) Tuning','Fontsize',16);
         plot([1,LO],[0,0],'k-');   
      end %End motion Tuning

      input 'check'
      close all
    end
      
end
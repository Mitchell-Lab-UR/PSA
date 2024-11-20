function PlotSimpleCriteria2(Info,Unit)
%******* 
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          Unit - spike cluster to run computation on

  %****** if Exp is [], then just perform the plot and return
  hh = figure('position',[100 100 800 400]); 
  HA = subplot('position',[0.10 0.15 0.35 0.7]);
  HB = subplot('position',[0.60 0.15 0.35 0.7]);
  %*******
  plot_stimpsth(Info,HA);
  plot_tuning(Info,HB);  
  
return;

%******* below is the plot of stimulus locked PSTH
function plot_stimpsth(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   
   tt = Info.StimTT;
   uu = Info.StimUU;
   su = Info.StimSU;
   basemu = Info.BaseMu;
   viscrit = Info.VisCrit;
   %********* now plot results
   maxo = max( uu + (2*su));
   %**** change to fill plot here
   h2 = plot(tt,uu,'k-',tt,(uu+(2*su)),'k-',tt,(uu-(2*su)),'k-'); hold on;
   %*** plot baseline rate line
   zz = find(tt < 0);  % get baseline firing as time before 0 ms
   if ~isnan(basemu)
      vmu = basemu;
      plot([Info.StimWin(1),Info.StimWin(2)],[vmu,vmu],'b--');
   end
   axis tight;
   V = axis;
   axis([-0.05 0.15 0 V(4)*1.05]);
   xlabel('Time (secs)');
   ylabel('Firing Rate (sp/s)');
   %******************
   
return;


%******* below is the plot of orientation tuning
function plot_tuning(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   OriTune = Info.OriTune;
   spikeplot.TuningCurvePlot(OriTune,[0,0,0.6],'.',0,0);
   %****
   axis tight;
   V = axis;
   axis([0 340 (V(3)*0.8) (V(4)*1.05)]);
   xlabel('Direction');
   ylabel('Rate');
   title(sprintf('DSI = %4.2f  conf( %4.2f to %4.2f )',OriTune.DSI,...
        (OriTune.DSI-(2*OriTune.DSI_SEM)),(OriTune.DSI+(2*OriTune.DSI_SEM))));    
  
return;
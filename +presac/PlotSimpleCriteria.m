function PlotSimpleCriteria(Info,Unit)
%******* 
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          Unit - spike cluster to run computation on

  %****** if Exp is [], then just perform the plot and return
  hh = figure('position',[100 100 800 400]); 
  HA = subplot('position',[0.20 0.17 0.30 0.7]);
  HB = subplot('position',[0.60 0.10 0.35 0.7]);
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
   zz = find( (tt>-0.01) & (tt<0.15));
   tt = tt(zz);
   uu = uu(zz);
   su = su(zz);
   %*******
   basemu = Info.BaseMu;
   viscrit = Info.VisCrit;
   %********* now plot results
   maxo = max( uu + (2*su));
   %**** change to fill plot here
   tta = [tt fliplr(tt)];
   uua = [(uu+(2*su)) fliplr(uu-(2*su))];
   fill(tta,uua,[0,0,0],'FaceAlpha',0.2,'Linestyle','none'); hold on;
   plot(tt,uu,'k-','Linewidth',2);
   %*** plot baseline rate line
   axis tight;
   V = axis;
   maxo = V(4)*1.05;
   axis([-0.02 0.12 0 maxo]);
   axis off;
   %*********
   plot([-0.01,0.12],[basemu,basemu],'k--');
   plot([0,0],[(0.02*maxo),(0.8*maxo)],'k-');
   text(-0.025,(0.9*maxo),'Onset','Fontsize',18);
   plot([0,0.05],[0,0],'k-','LineWidth',2);
   text(0.01,-5,'50 ms','Fontsize',20);
   plot([-0.02,-0.02],[0,25],'k-','LineWidth',2);
   text(-0.04,10,'25(sp/s)','Fontsize',20,'Rotation',90);
   %******************
   
return;


%******* below is the plot of orientation tuning
function plot_tuning(Info,H)
   %******* This plots 1) firing rate as a function of orientation with error bars
   
   subplot(H);
   OriTune = Info.OriTune;
   spikeplot.TuningCurvePlot2(OriTune,[0,0,0],'.',0,0);
   %****
   axis tight;
   V = axis;
   maxo = V(4)*1.05;
   axis([-30 340 -5 maxo]);
   if (1)
      plot([-15,-15],[0,20],'k-','LineWidth',2);
      text(-50,0,'20(sp/s)','Fontsize',20,'Rotation',90);
      plot([0,340],[0,0],'k-','Linewidth',2);
      plot([0,0],[0,-2],'k-','Linewidth',2);
      plot([340,340],[0,-2],'k-','Linewidth',2);
      text(-10,-5,'0','Fontsize',14);
      text(330,-5,'340','Fontsize',14);
      text(40,-7,'Direction(degrees)','Fontsize',18);
   end 
   axis off;
   text(60,(0.2*maxo),sprintf('DSI=%4.2f(+/- %4.2f)',OriTune.DSI,...
            (2*OriTune.DSI_SEM)),'Fontsize',18);
  
return;
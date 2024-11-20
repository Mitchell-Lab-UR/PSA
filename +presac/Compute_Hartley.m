function Info = Compute_Hartley(Info,Exp,Unit,H)
%******* 
%****** Based on: function Info = Compute_DSI(Info,Exp,Unit,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          Unit - spike cluster to run computation on
%***          H - handle to plot window if provided, or if 1 create, else []
%***
%*** Outputs: Info - it computes the orientation tuning and spatial frequency of all trials

Info.Orientation=[];
Info.SpatialFrequency=[];
  %****** if Exp is [], then just perform the plot and return
  if isempty(Exp)
    if isfloat(H) && (H > 0)
       hh = figure('position',[100 100 1200 400]); 
       HA = subplot('position',[0.0752 0.15 0.25 0.7]);
       HB = subplot('position',[0.4 0.15 0.25 0.7]);
       HC = subplot('position',[0.725 0.15 0.25 0.7]);
    else
        if isempty(H)
          return;
        else
          if (length(H) ~= 2)
              disp('Error with subplot pass to Compute_Hartley, requires two subplots');
              return;
          else
              HA = H{1};
              HB = H{2};
              HC = H{3};
          end
        end
    end
    %*******
%     plot_compute_hartley(Info,HA);  % temporal kernel
%     plot_orientation(Info,HB);      % Orientation
%     plot_spatialfrequency(Info,HC); % Spatial Frequency

% looks at temporal kernel to gratings, plus spat freq x orientation
SPClust = 1;
[StimX,StimY] = Forage.StimMatrix_ForageGratingKernel(Exp,SPClust);  %step 1
ForageGratingKernel_Ori;     % step 2, RF

    if (size(ORI)) <= 2
        Info.Orientation=[ORI];
    else
        Info.Orientation= NaN;
    end
    
    if (size(SpatFreq)) <= 1
        
        Info.SpatialFrequency=SpatFreq;
    else
        Info.SpatialFrequency= NaN;
    end
    
    %********
    if isfloat(H)  && (H == 2)
       z = getframe(hh);  % current figure
       uname = [Info.pathplot,filesep,'Hartley_',Info.tagname,'.png'];
       disp(sprintf('Storing image of Hartley graph at %s',uname)); 
       imwrite(z.cdata,uname); % store image for later review
       close(hh);
    end   
    return;
  end
  

%******* below is the plot function in the same file
function plot_compute_hartley(Info,H)
   subplot(H);
   plot(tXX,tempker,'k-'); hold on;
   plot(tXX,zeros(size(tempker)),'k--');
return;

%******* below is the plot of orientation tuning
function plot_orientation(Info,H)
   %******* This plots orientation average
    uu = reshape(stcounts,[NOri NSpf]);
    su = reshape(stsem,[NOri NSpf]);   
    otune = mean(uu');
    otune = [otune otune];
    sotune = mean(su');
    sotune = [sotune sotune];
    xo = [SpatOris (SpatOris+180)];
    plot(xo,otune,'bo'); hold on;
    plot(xo,otune + (2*sotune),'b--'); hold on;
    plot(xo,otune - (2*sotune),'b--'); hold on;
    plot(xo,ymean*ones(size(xo)),'k-');
    axis tight;
    V = axis;
    mag = 0.1*(V(4)-V(3));
    axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
    xlabel('Orientation (degs)');
    ylabel('Rate (sp/s)');
   
return;


%******* below is the plot of SpatFreq
function plot_spatialfrequency(Info,H)
   %******* This plots spatail frequency tuning
    uu = reshape(stcounts,[NOri NSpf]);
    su = reshape(stsem,[NOri NSpf]); 
    stune = mean(uu);
    sstune = mean(su);
    xo = SpatFrqs; 
    semilogx(xo,stune,'bo'); hold on;
    semilogx(xo,stune + (2*sstune),'b--'); hold on;
    semilogx(xo,stune - (2*sstune),'b--'); hold on;
    semilogx(xo,ymean*ones(size(xo)),'k-');
    axis tight;
    V = axis;
    mag = 0.1*(V(4)-V(3));
    axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
    xlabel('SpatFreq (cyc/deg)');
    ylabel('Rate (sp/s)');
    set(gca,'XTick',MyTick,'XTickLabel',MyTickLabel);
   
   %***************
   
return;
   
function Info = Compute_SplitPFR(Info,Exp,H)
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
       %****** setup a multi-panel plot (four panels) 
       hh = figure('position',[400 200 1400 400]); 
       HA = subplot('position',[0.075 0.15 0.18 0.7]);
       HB = subplot('position',[0.3 0.15 0.18 0.7]);
       HC = subplot('position',[0.55 0.15 0.18 0.7]);
       HD = subplot('position',[0.8 0.15 0.18 0.7]);
       
       %********************
       plot_PFR_unity_Best(Info,HA);   % plot true vs estimate direction
       plot_PFR_unity_Worst(Info,HB);   % plot true vs estimate direction
       plot_PFR_angles_Best(Info,HC);  % plot of angle distribution
       plot_PFR_angles_Worst(Info,HD);  % plot of angle distribution
       if (H == 2)  % store result to png and close
          z = getframe(hh);  % current figure
          uname = [Info.pathplot,filesep,'PFR_',Info.tagname,'.png'];
          disp(sprintf('Storing image of PFR graph at %s',uname)); 
          imwrite(z.cdata,uname); % store image for later review
          close(hh);
       end
    else 
       if ~isempty(H)  % H should be a sub-plot handle otherwise
          plot_PFR_unity_Best(Info,H{1});   % plot true vs estimate direction
          plot_PFR_unity_Worst(Info,H{2});   % plot true vs estimate direction
          plot_PFR_angles_Best(Info,H{3});  % plot of angle distribution
          plot_PFR_angles_Worst(Info,H{4});  % plot of angle distribution
       end
    end
    return;
  end

  NT = Info.NumTargs;
  
  %****** compute median angular error for target
  %****** and for the non-targets per session (for best and worst trials)
  dots = Info.pfrs{1}(:,1); 
  cros = Info.pfrs{1}(:,2); 
  tglength =  Info.OpenLoopSpeed * (Info.OpenLoopEnd-Info.OpenLoopStart);
  mago = abs( dots + i * cros);
  zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
  angos = angle( dots(zz) +  i * cros(zz) ) * (180/pi);  % go from radians to degs
  
   % find where to spilt the data 
   med_angos = median(abs(angos));
   
   yy = find( abs(angos) < med_angos);  % take Best PFR
   %********
   best_angos = angle( dots(yy) +  i * cros(yy) ) * (180/pi);  % go from radians to degs
  Info.PFR_TargErr_Best = median(abs(best_angos));
  
  xx = find( abs(angos) >= med_angos);  % take Worst PFR
   %********
   worst_angos = angle( dots(xx) +  i * cros(xx) ) * (180/pi);  % go from radians to degs
  Info.PFR_TargErr_Worst = median(abs(worst_angos));
 
  
  %*** resultant vector (Best PFR)
  udot = mean(dots(yy));
  ucro = mean(cros(yy));
  Info.PFR_TargResV_Best = [udot,ucro];
  %**** compute histograms of PFR error
  vx = -180:10:180;
  vy = hist(best_angos,vx);
  vy = vy / sum(vy);   % normalize
  %*******
  Info.PFR_TargHistX_Best = vx;
  Info.PFR_TargHist_Best = vy;
  
  %*** resultant vector (Worst PFR)
  udot = mean(dots(xx));
  ucro = mean(cros(xx));
  Info.PFR_TargResV_Worst = [udot,ucro];
  %**** compute histograms of PFR error
  vx = -180:10:180;
  vy = hist(worst_angos,vx);
  vy = vy / sum(vy);   % normalize
  %*******
  Info.PFR_TargHistX_Worst = vx;
  Info.PFR_TargHist_Worst = vy;
  %*******
  
  
    %******* Now Calculate PRF for distractors 
  if (NT > 2)
    dots = [ Info.pfrs{2}(yy,1); Info.pfrs{3}(yy,1)];
    cros = [ Info.pfrs{2}(yy,2); Info.pfrs{3}(yy,2)];
  else
    dots = [ Info.pfrs{2}(yy,1)];
    cros = [ Info.pfrs{2}(yy,2)]; 
  end
  mago = abs( dots + i * cros);
  zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
  angos = angle( dots(zz) +  i * cros(zz) ) * (180/pi);  % go from radians to degs
  Info.PFR_DistErr_Best = median(abs(angos));
  vy = hist(angos,vx);
  vy = vy / sum(vy);   % normalize
  
  %*** resultant vector
  udot = mean(dots(zz));
  ucro = mean(cros(zz));
  Info.PFR_DistResV_Best = [udot,ucro];
  %*******
  Info.PFR_DistHist_Best = vy;
   
  if (NT > 2)
    dots = [ Info.pfrs{2}(xx,1); Info.pfrs{3}(xx,1)];
    cros = [ Info.pfrs{2}(xx,2); Info.pfrs{3}(xx,2)];
  else
    dots = [ Info.pfrs{2}(xx,1)];
    cros = [ Info.pfrs{2}(xx,2)];    
  end
  mago = abs( dots + i * cros);
  zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
  angos = angle( dots(zz) +  i * cros(zz) ) * (180/pi);  % go from radians to degs
  Info.PFR_DistErr_Worst = median(abs(angos));
  vy = hist(angos,vx);
  vy = vy / sum(vy);   % normalize
  %*** resultant vector
  udot = mean(dots(zz));
  ucro = mean(cros(zz));
  Info.PFR_DistResV_Worst = [udot,ucro];
  %*******
  Info.PFR_DistHist_Worst = vy;
  
  
  
  
  
  
  
  
  
  
   
return;


%***** plot the true motion versus estimated direction on a unity line          
function plot_PFR_unity_Best(Info,H)
   % pfr directions
   dots = Info.pfrs{1}(:,1); 
   cros = Info.pfrs{1}(:,2); 
   tangs = Info.pfrs{1}(:,4); % actual motion of target apertures
   %*** throw out the outliers (huge PFR due to 2nd sac)
   tglength =  Info.OpenLoopSpeed * (Info.OpenLoopEnd-Info.OpenLoopStart);
   mago = abs( dots + i * cros);
   zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
   %********
   angos = angle( dots(zz) +  i * cros(zz) ) * (180/pi);  % go from radians to degs
   tangs = tangs(zz); % match vectors
   %*********
   % find where to spilt the data 
   med_angos = median(abs(angos));
   
   yy = find( abs(angos) < med_angos);  % take subset of total with reasonable PFR
   %********
   angos = angle( dots(yy) +  i * cros(yy) ) * (180/pi);  % go from radians to degs
   tangs = tangs(yy); % match vectors
   
   subplot(H); 
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
   plot(tangs,angos,'g.'); hold on;
   plot([0,360],[0,360],'k-'); 
   axis([0 360 0 360]);
   xlabel('True motion direction');
   ylabel('Best PFR motion direction');
   title(Info.tagname);
   %****
   
return;

function plot_PFR_unity_Worst(Info,H)
   % pfr directions
   dots = Info.pfrs{1}(:,1); 
   cros = Info.pfrs{1}(:,2); 
   tangs = Info.pfrs{1}(:,4); % actual motion of target apertures
   %*** throw out the outliers (huge PFR due to 2nd sac)
   tglength =  Info.OpenLoopSpeed * (Info.OpenLoopEnd-Info.OpenLoopStart);
   mago = abs( dots + i * cros);
   zz = find( mago < (2*tglength));  % take subset of total with reasonable PFR
   %********
   angos = angle( dots(zz) +  i * cros(zz) ) * (180/pi);  % go from radians to degs
   tangs = tangs(zz); % match vectors
   %*********
   % find where to spilt the data 
   med_angos = median(abs(angos));
   
   yy = find( abs(angos) >= med_angos);  % take subset of total with reasonable PFR
   %********
   angos = angle( dots(yy) +  i * cros(yy) ) * (180/pi);  % go from radians to degs
   tangs = tangs(yy); % match vectors
   
   subplot(H); 
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
   plot(tangs,angos,'r.'); hold on;
   plot([0,360],[0,360],'k-'); 
   axis([0 360 0 360]);
   xlabel('True motion direction');
   ylabel('Worst PFR motion direction');
   title(Info.tagname);
   %****
   
return;

function plot_PFR_angles_Best(Info,H)  % plot of angle distribution
          
   subplot(H);
   %**********
   plot(Info.PFR_TargHistX_Best,Info.PFR_TargHist_Best,'g.-'); hold on;
   plot(Info.PFR_TargHistX_Best,Info.PFR_DistHist_Best,'b.-'); hold on;
   %***********
   xlabel('Angle (degs)');
   ylabel('Counts');
   mederr = Info.PFR_TargErr_Best;
   title(sprintf('(Best) Med Targ Err: %5.1f',mederr));
   %**********
   
return;

function plot_PFR_angles_Worst(Info,H)  % plot of angle distribution
          
   subplot(H);
   %**********
   plot(Info.PFR_TargHistX_Worst,Info.PFR_TargHist_Worst,'r.-'); hold on;
   plot(Info.PFR_TargHistX_Worst,Info.PFR_DistHist_Worst,'b.-'); hold on;
   %***********
   xlabel('Angle (degs)');
   ylabel('Counts');
   mederr = Info.PFR_TargErr_Worst;
   title(sprintf('(Worst) Med Targ Err: %5.1f',mederr));
   %**********
   
return;




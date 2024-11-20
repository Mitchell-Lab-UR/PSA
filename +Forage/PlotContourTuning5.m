function [CT,mloc] = PlotContourTuning5(rfinfo,PlotIt,BaseTag)
   % based on function PlotStackedRF_Summary(rfinfo,PlotIt,BaseTag)  
   %**** loops through RFs and plots their contours
   % rfinfo - struct with RF info for all units
   % PlotIt - show plot or not
   % BaseTag - name of source file to attach in plot
   %
   % Outputs:  CT - contour of x,y points, and other information
   %           mloc - median centroid of contour points
   
   CT = [];
   mloc = [];
   if ~isempty( rfinfo )
       %**********
       if (sum(sum(rfinfo.sigcounts)) < 5)
           return;  % require a significant spatial RF
       end
       %*********     
       mino = mean(mean(rfinfo.meanrf)); %RFPlotMino; %min(min(rfinfo{k}.mcounts));
       maxo = max(max(rfinfo.meanrf)); %rfinfo{k}.RFPlotMaxo; %max(max(rfinfo{k}.mcounts));
       %*** plot a contour line here
       watermarks = [0.35,0.5,0.65,0.8,0.95];
       ID_mark = 2;  % which watermark to ID basic RF contour
       ID_top = 4;  % top contour for peaks
       NW = length(watermarks);
       cx = cell(1,NW);
       cy = cell(1,NW);
       cpoly = cell(1,NW);
       for ck = 1:NW 
          v = mino + watermarks(ck)*(maxo-mino);  % half-height of RF     
          [cx{ck},cy{ck},cpoly{ck}] = compute_contour_rf(rfinfo.meanrf,rfinfo.Zx,rfinfo.Zy,v);
       end
       %***** use ID_mark to find centroid and define core RF area
       cpoly{ID_mark} = compute_dominant_contour(cpoly{ID_mark},[]);  % find single dominant
       [xx,yy] = centroid(cpoly{ID_mark});
       mloc = [xx,yy];
       %****** CT is still empty:  before you load anything into CT,
       %**** determine if the "marked" dominnant boundary is closed ...
       %**** and if not (part is off the screen), then return CT=[]
       verts = cpoly{ID_mark}.Vertices;
       verts = [verts ; verts(1,:)];  % bound circle
       dverts = diff(verts);
       for k = 1:length(dverts)
           ndverts(k) = norm(dverts(k,:));
       end
       if (length(ndverts) < 10) || (max(ndverts) > (5*median(ndverts)))
         disp('Identified 5 times gap in main contour, not closed?? Skipping');
         CT = [];
         return;
       end
       %****** store basic results
       CT.cx = cx{ID_mark};
       CT.cy = cy{ID_mark};
       CT.cx2 = cx{ID_top};
       CT.cy2 = cy{ID_top};
       CT.mloc{1} = mloc;
       CT.watermarks = watermarks;
       CT.ID_mark = ID_mark;
       CT.ID_top = ID_top;
       %****** Now loop on others and use intersect to find contours
       %****** to the current contour, throwing away others, and if more
       %****** than one intersecting, taking two at maximum
       for ck = 1:NW
           if (ck == ID_mark)
               continue;
           else
               if (ck < ID_mark)
                 cpoly{ck} = compute_contact_contour(cpoly{ck},cpoly{ID_mark});  % find intersects with single dominant
               else
                 cpoly{ck} = compute_dominant_contour(cpoly{ck},cpoly{ID_mark});  % find intersects with single dominant
               end           
           end
       end
       CT.cpoly = cpoly;
     
       %**** simplify shapes to encode them
       NLEN = 50;
       NCT = cell(1,NW);
       RCT = cell(1,NW);
       [NCT{ID_mark},RCT{ID_mark}] = normed_polyshape(cpoly{ID_mark},mloc,NLEN);
       for ck = 1:NW
           if (ck == ID_mark)
               continue;
           else   
               if (cpoly{ck}.NumRegions > 1)
                  %** store the best two regions in proportion to their size
                  pg = sortregions(cpoly{ck},'numsides','descend');
                  R = regions(pg);
                  cpolyA = R(1);
                  cpolyB = R(2);
                  NA = length(cpolyA.Vertices);
                  NB = length(cpolyB.Vertices);
                  NLENA = floor( (NLEN-1)*(NA/(NA+NB)) );
                  NLENB = floor( (NLEN-1)*(NB/(NA+NB)) );
                  [NCT2A,RCT2A] = normed_polyshape(cpolyA,mloc,NLENA); 
                  [NCT2B,RCT2B] = normed_polyshape(cpolyB,mloc,NLENB); 
                  GAP = (NLEN-NLENA-NLENB);
                  NCT{ck} = [NCT2A ; nan(GAP,2) ; NCT2B];
                  RCT{ck} = [RCT2A ; nan(GAP,2) ; RCT2B];
                  %********* dual region detectable by NaN inside it
                  if (ck == ID_top)
                     [xx,yy] = centroid(cpolyA);
                     CT.mloc{2} = [xx,yy];
                     [xx,yy] = centroid(cpolyB);
                     CT.mloc{3} = [xx,yy];
                  end
               else
                  [NCT{ck},RCT{ck}] = normed_polyshape(cpoly{ck},mloc,NLEN); 
                  if (ck == ID_top)
                     [xx,yy] = centroid(cpoly{ck});
                     CT.mloc{2} = [xx,yy];
                     CT.mloc{3} = [NaN,NaN];
                  end
               end 
           end
       end
       %******** make total RF of many layers for plotting
       CT.NCT = NCT;
       CT.RCT = RCT;
       CT.Nv = [];
       CT.Rv = []; 
       for ck = 1:NW
          CT.Nv = [ CT.Nv; [NaN,NaN]; NCT{ck}];
          CT.Rv = [ CT.Rv; [NaN,NaN]; RCT{ck}];
       end
       
       
       %******* get motion tuning
       ND = rfinfo.ND;
       xx = (0:(ND+6))*(360/ND);
       yy = [1:ND,1:(size(xx,2)-ND)];
       CT.mo_x = xx;
       CT.mo_u = rfinfo.mou(yy);
       CT.mo_std = rfinfo.mostd(yy);
       CT.sacbase = rfinfo.sacbase;
       CT.postcurve = rfinfo.postcurve(yy);
       CT.postcurve_std = rfinfo.postcurve_std(yy);
       CT.precurve = rfinfo.precurve(yy);
       CT.precurve_std = rfinfo.precurve_std(yy);
       CT.durcurve = rfinfo.durcurve(yy);
       CT.durcurve_std = rfinfo.durcurve_std(yy);
       %**************
       SACN = rfinfo.SACN;
       CT.sac_x = rfinfo.tSXX;
       CT.sac_u = rfinfo.mvec(SACN+1,:);
       CT.sac_std = rfinfo.mvec2(SACN+1,:);
       sac_x = CT.sac_x;
       sac_u = CT.sac_u;
       sac_std = CT.sac_std;
       sbase = mean( sac_u( find( sac_x < 0 ) ) );  % baseline rate
       sac_z = (CT.sac_u - sbase) ./ sac_std; 
       sac_z = smooth(sac_z,5)';
       sac_z = sbase + (sac_z * mean(sac_std));
       CT.sac_z = sac_z;
       %********
       
       %*** DEFINE a CT.stats as a return of needed stats
       %** Columns:  RF size, RF height, RF width, RF elongation, DSI, ...
       %****            SacMod, SacSuppression, SacRebound,...
       %****            SacMinTime,SacMaxTime
       %*********
       HalfMark = find( watermarks == 0.5);
       nct = CT.NCT{HalfMark};
       rct = CT.RCT{HalfMark};
       RF_size = nanmean(rct(:,2));  % mean radius
       RF_height = range(nct(:,2));  % max to min on y-axis 
       RF_width = range(nct(:,1));  % max to min on y-axis 
       RF_elongation = RF_height / RF_width;
       %*** compute DSI
       mo_x = CT.mo_x;
       mo_u = CT.mo_u;
       zz = find( (mo_x >= 0) & (mo_x < 360));
       wadd = complex(0,0);
       wsum = 0;
       for ak = 1:length(zz)
           ango = mo_x(zz(ak)) * (pi/180);
           rate = mo_u(zz(ak));
           wadd = wadd + complex( rate*cos(ango), rate*sin(ango));
           wsum = wsum + rate;
       end
       if (wsum)
           DSI = abs( wadd/wsum );
       else
           DSI = NaN;
       end
       %****** saccade modulation stats
       sac_x = CT.sac_x;
       sac_u = CT.sac_u;
       sac_std = CT.sac_std;
       sac_z = CT.sac_z;
       %*****
       %******
       zz = find( (sac_x > 0) & (sac_x < 100));
       zzmin = find( sac_z(zz) == min(sac_z(zz)));
       SacMinTime = sac_x(zz(zzmin(1)));
       zz = find( (sac_x > 0) & (sac_x < 200));
       zzmax = find( sac_z(zz) == max(sac_z(zz)));
       SacMaxTime = sac_x(zz(zzmax(1)));
       %***
       Fudge = 15;
       zzminset = find( (sac_x >= (SacMinTime-Fudge)) & (sac_x < (SacMinTime+Fudge)) )
       zzmaxset = find( (sac_x >= (SacMaxTime-Fudge)) & (sac_x < (SacMaxTime+Fudge)) )
       smax = mean( sac_u(zzmaxset) );
       smin = mean( sac_u(zzminset) );
       SacMod = (smax-smin)/(smax+smin);
       SacSuppression = (sbase-smin)/(sbase+smin);
       SacRebound = (smax-sbase)/(smax+sbase);
       %*********
       CT.stats = [RF_size,RF_height,RF_width,RF_elongation,DSI,...
                      SacMod,SacSuppression,SacRebound,SacMinTime,SacMaxTime];
       CT.stats           
       %***********
       
       if (PlotIt)
           gcol = [0,0.6,0.298];
           rad = norm(mloc);
           Limo = (3*rad);
           %**********
           hf = figure;
           set(hf,'position',[800 250 800 800]);
           %*******
           subplot('position',[0.1 0.55 0.30 0.30]);
           plot(CT.cx,CT.cy,'k.-'); hold on;
           h = plot(CT.cx2,CT.cy2,'k.-'); hold on;
           set(h,'Color',[0.5,0.5,0.5]);
           plot(CT.mloc{1}(1),CT.mloc{1}(2),'k+');  
           plot(CT.mloc{2}(1),CT.mloc{2}(2),'ks');  
           plot(CT.mloc{3}(1),CT.mloc{3}(2),'ko');  
           plot([-Limo,Limo],[0,0],'k:');
           plot([0,0],[-Limo,Limo],'k:');
           axis([-Limo Limo -Limo Limo]);
           xlabel('Horizontal (degs)');
           ylabel('Vertical (degs)');
           title('Contour Original Space');

           %******* show stacked in grey scale (pretty figure)
           rmaxo = max( max(max(CT.Rv(:,2))), 1.5);
           rmaxo = max(max(CT.Rv(:,2)));
           if (rmaxo < 1.5)
               rmaxo = 1.5;
           end
           pax = polaraxes('position',[0.6 0.55 0.30 0.30],'RLim',[0 rmaxo]);
           for ck = 1:length(CT.RCT)
               val = (length(CT.RCT)-ck)/length(CT.RCT);
               h = polarplot(pax,CT.RCT{ck}(:,1)',CT.RCT{ck}(:,2)','k-'); hold on;
               set(h,'Color',[val,val,val]);
               set(h,'LineWidth',(1+(2*(1-val))));
           end
           theta = 0:0.1:(2*pi);
           hg = polarplot(pax,theta,ones(size(theta)),'r-');
           set(hg,'Color',gcol);
           pax.RLim = [0 rmaxo];     
           title('Normalized RF by Contours'); 
           %**********
 
           %******* Show motion tuning, and saccade modulation by direction
           subplot('position',[0.1 0.15 0.30 0.30]);
           h = plot(CT.mo_x,CT.mo_u,'k-'); hold on;
           set(h,'LineWidth',2);
           plot(CT.mo_x,CT.mo_u+(2*CT.mo_std),'k-'); hold on;
           plot(CT.mo_x,CT.mo_u-(2*CT.mo_std),'k-'); hold on;
           if(0) % plot sac direction or not
           plot([CT.mo_x(1),CT.mo_x(end)],[CT.sacbase,CT.sacbase],'k-');
           %****
           h = plot(CT.mo_x,CT.postcurve,'b--'); hold on;
           set(h,'LineWidth',2);
           % plot(CT.mo_x,CT.postcurve+(2*CT.postcurve_std),'b-'); hold on;
           % plot(CT.mo_x,CT.postcurve-(2*CT.postcurve_std),'b-'); hold on;
           %****
           h = plot(CT.mo_x,CT.precurve,'r--'); hold on;
           set(h,'LineWidth',2);
           % plot(CT.mo_x,CT.precurve+(2*CT.precurve_std),'r-'); hold on;
           % plot(CT.mo_x,CT.precurve-(2*CT.precurve_std),'r-'); hold on;
           %****
           h = plot(CT.mo_x,CT.durcurve,'m--'); hold on;
           set(h,'LineWidth',2);
           % plot(CT.mo_x,CT.durcurve+(2*CT.durcurve_std),'m-'); hold on;
           % plot(CT.mo_x,CT.durcurve-(2*CT.durcurve_std),'m-'); hold on;
           %*******  
           end
           xlabel('Direction (degs)');
           ylabel('Rate (hz)');
           title('Motion and Saccade(pre,dur,post) Tuning');
           
           %******** Show temporal modulation by saccade
           subplot('position',[0.6 0.15 0.30 0.30]);
           h = plot(CT.sac_x,CT.sac_u,'k-'); hold on;
           set(h,'LineWidth',2);
           plot(CT.sac_x,CT.sac_u+(2*CT.sac_std),'k-'); hold on;
           plot(CT.sac_x,CT.sac_u-(2*CT.sac_std),'k-'); hold on;
           plot(CT.sac_x,CT.sac_z,'r-'); hold on;
           plot([CT.sac_x(1),CT.sac_x(end)],[CT.sacbase,CT.sacbase],'k-');
           xlabel('Time (ms)');
           ylabel('Rate (hz)');
           title('Saccade Modulation');
           %*********** 
       end            
   end
   
return;


function  [NCT,RCT] = normed_polyshape(cpoly2,mloc,NLEN)

       %******* ROTATE CENTER OF SHAPE TO ALIGN ON UP AXIS
       ango = angle(complex(mloc(1),mloc(2)));
       dango = ango-(pi/2);   % rotate to 180 degs
       ROTO = [[cos(dango) sin(dango)]; [-sin(dango) cos(dango)]];
       rr = (ROTO * cpoly2.Vertices')';
       scal = 1/norm(mloc);
       rr = rr * scal;   % renorm to fix at 1 deg ecc (normalized RF)
       
       %******* now go around the contour in equal spaced steps
       RN = length(rr);
       %****** first, find 0 degrees, and how to step clockwise
       GOPTS = [];
       for theta = [0,30]  % find start and direction to move
          ango = theta * (180/pi);
          v = [cos(ango),sin(ango)];
          %****** find the closest contour point to that angle
          best = -Inf;
          bestk = NaN;
          bestrad = 0;
          for k = 1:RN
              dv = rr(k,:) - [0,1];
              rad = norm(dv);
              dv = dv / rad;
              doto = sum( dv .* v);
              if (doto > best)
                  best = doto;
                  bestk = k;
                  bestrad = rad;
              end
          end
          GOPTS = [GOPTS ; bestk];
       end  
       startk = GOPTS(1);
       deltak = sign( GOPTS(2) - GOPTS(1));
       %**** compute distance around contour on those steps
       TOTDIST = 0;
       cmdist = 0;
       cmk = startk;
       k = startk;
       tally = 0;
       while (tally < RN)
           %******* find next step
           newk = k + deltak;
           if (newk <= 0)
                 newk = RN;
           end
           if (newk > RN)
                 newk = 1; 
           end
           dist = norm(rr(newk,:)-rr(k,:));
           TOTDIST = TOTDIST + dist;
           cmdist = [cmdist ; TOTDIST];
           cmk = [cmk ; newk];
           k = newk;
           %*******
           tally = tally + abs(deltak);
       end    
       %********* interpolate the points of equal spacing around contour
       STEPDIST = TOTDIST/(NLEN-2);
       NCT = [];
       RCT = [];
       for k = 1:(NLEN-1)
           val = (k-1)*STEPDIST;
           zza = find(cmdist <= val);
           if isempty(zza)
               ka = cmk(1);
               wkb = 0;
           else
               ka = cmk(zza(end));
               wkb = (val-cmdist(ka));
           end
           zzb = find(cmdist > val);
           if isempty(zzb)
               kb = cmk(end);
               wka = 0;
           else
               kb = cmk(zzb(1));
               wka = (val-cmdist(ka));
           end
           loco = (wka*rr(ka,:) + wkb*rr(kb,:))/(wka+wkb);  % linear interp
           %*******
           NCT = [NCT ; loco];
           mdiff = loco-[0,1];  % subtrack normed center
           ango = angle(complex(mdiff(1),mdiff(2)));
           brad = norm(mdiff);
           RCT = [RCT ; [ango brad]];
           %**************
       end
       NCT = [NCT ; NCT(1,:)];  % force self loop
       RCT = [RCT ; RCT(1,:)];  % force self loop
       
 return;

 
 function [cx,cy,cpoly] = compute_contour_rf(mrf,Zx,Zy,v,vpoly)
       
       %******* upsample the RF before getting contours? *****
       [X,Y] = meshgrid(1:size(mrf,1),1:size(mrf,2));
       [Xq,Yq] = meshgrid(1:0.1:size(mrf,1),1:0.1:size(mrf,2));
       meanrf = interp2(X,Y,mrf,Xq,Yq,'spline');
       %***********
       c = contourc(meanrf,[v v]);
       cx = c(1,2:end);
       cy = c(2,2:end);
       %******* only accept smooth contours
       for ik = 2:length(cx)
           dist = norm(cx(ik)-cx(ik-1),cy(ik)-cy(ik-1));
           if (dist > 4)
               cx((ik-1):ik) = NaN;
               cy((ik-1):ik) = NaN;
           end
       end
       Nx = size(meanrf,1);
       Ny = size(meanrf,2);
       cx = ((cx-1)/(Nx-1))*range(Zx) + min(Zx);
       cy = ((cy-1)/(Ny-1))*range(Zy) + min(Zy);
       %**********
       xmin = min(Zx);
       xmax = max(Zx);
       ymin = min(Zy);
       ymax = max(Zy);
       %*******
       cy = -cy;  % need to invert the y
 
       %******* filter down to dominant set of vertices
       cpoly = polyshape(cx',cy');
      
 return;
 
 function cpoly = compute_dominant_contour(cpoly,vpoly)
       %******* filter down to dominant set of vertices
       %**** vpoly should be a single regions, get cpoly inside it
       if isempty(vpoly)
         pg = sortregions(cpoly,'numsides','descend');
         R = regions(pg);
         cpoly = R(1);
       else
          cpoly = intersect(vpoly,cpoly);  % take anything inside vpoly
       end        
 return;
         
 function zpoly = compute_contact_contour(cpoly,vpoly)
       % find the subregion in cpoly that overlaps vpoly most      
       zpoly = [];
       pg = sortregions(cpoly,'numsides','descend');
       R = regions(pg);
       bestk = NaN;
       bestaa = 0;
       for k = 1:length(R)
          newpoly = intersect(R(k),vpoly);
          aa = area(newpoly);
          if (aa > bestaa)
              bestaa = aa;
              bestk = k;
          end
       end
       if ~isnan(bestk)
           zpoly = R(bestk);
       end
 return;
       
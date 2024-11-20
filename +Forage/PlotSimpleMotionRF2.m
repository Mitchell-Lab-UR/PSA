function [cx,cy] = PlotSimpleMotionRF2(rfinfo)
% function PlotForageSpatialKernel(rfinfo,ftag)
%   shows the results of RF analysis as plots, show file tag for unit name

%******* take all fields of rfinfo and make part of the environment
%******** download variables stored in info
fields = fieldnames(rfinfo);
for k = 1:size(fields,1)
      str = [fields{k} ' = rfinfo.' fields{k} ';'];
      eval(str);
end
%*******************************

    hf = figure;
    set(hf,'position',[900 100 600 600]);  
    
    %**** plot RF across time frames, and show significant points
    istart = 2;
    iend = 4;
    IKN = (KN-iend);
    %*********
    dx = (0.6/(IKN-istart));
    dxs = (0.6/(IKN-istart+1));
    mino = RFPlotMino; %min(min(mcounts));
    maxo = RFPlotMaxo; %max(max(mcounts));
    for it = istart:IKN
       fx = 0.30 + ((it-istart)*dx);
       subplot('position',[0.15 (1-fx) dxs dxs]);
       ito = (2-DTA)+it;
       svec = flipud( reshape(squeeze(mcounts(:,ito)),Nx,Ny)' );
       svec = imgaussfilt(svec,1.0);
       hi = imagesc(Zx,Zy,svec,[mino maxo]); hold on;    
       % plot([-15,15],[0,0],'k-');
       % plot([0,0],[-15,15],'k-');
       axis off;
       % h = xlabel(sprintf('%4.1f ms',tXX(ito)),'Fontsize',12);   
       h = text(-Nx*1.5,Ny*0.7,sprintf('%2.0f',tXX(ito)),...
                 'Fontsize',16,'Rotation',90);
       if (it == IKN) % floor(((1*istart)+(5*IKN))/6) )
            h = text(-Nx*2.8,-Ny*0.5,'Delay from dot onset (ms)',...
                 'Fontsize',24,'Rotation',90);
       end
    end
    %***********
    subplot('position',[0.3 0.25 0.6 0.5]);
    smeanrf = imgaussfilt(meanrf,0.5);
    imagesc(Zx,Zy,smeanrf,[mino maxo]); hold on;
    plot([-15,15],[0,0],'k:','Linewidth',2);
    plot([0,0],[-15,15],'k:','Linewidth',2);
    plot([0,14.5],[14,14],'k-','Linewidth',3,'Color',[0.99,0.99,0.99]);
    ht = text(3,11,'15 dva','Fontsize',20,'Color',[0.99,0.99,0.99]);
    axis off;
    % title('Spatial RF','Fontsize',18);
    
    %******
    mino = mean(mean(smeanrf)); 
    maxo = max(max(smeanrf)); 
    %*** plot a contour line here
    v = mino + 0.5*(maxo-mino);  % half-height of RF
    c = contourc(smeanrf,[v v]);
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
    %******
    cx = ((cx-1)/(Nx-1))*range(Zx) + min(Zx);
    cy = ((cy-1)/(Ny-1))*range(Zy) + min(Zy);
    hh = plot(cx,cy,'k-'); hold on;
    %*****
    hc = colorbar;
    set(hc,'Fontsize',14);
    hc2 = text(Nx*1.6,Ny*0.6,'Firing Rate (sp/s)','Fontsize',20,'Rotation',90);
    text(-Nx*0.8,-Ny*1.2,'RF at Peak Delay','Fontsize',24);
    %************

return;

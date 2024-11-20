function [cx,cy] = PlotSimpleMotionRF(rfinfo)
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
    set(hf,'position',[150 50 1000 333]);  
    
    %**** plot RF across time frames, and show significant points
    istart = 2;
    iend = 4;
    IKN = (KN-iend);
    %*********
    dx = (0.45/(IKN-istart));
    dxs = (0.45/(IKN-istart+1));
    mino = RFPlotMino; %min(min(mcounts));
    maxo = RFPlotMaxo; %max(max(mcounts));
    for it = istart:IKN
       fx = 0.025 + ((it-istart)*dx);
       subplot('position',[fx 0.25 dxs (3.2*dxs)]);
       ito = (2-DTA)+it;
       svec = flipud( reshape(squeeze(mcounts(:,ito)),Nx,Ny)' );
       svec = imgaussfilt(svec,1.0);
       hi = imagesc(Zx,Zy,svec,[mino maxo]); hold on;    
       % plot([-15,15],[0,0],'k-');
       % plot([0,0],[-15,15],'k-');
       axis off;
       h = title(sprintf('%4.1f ms',tXX(ito)),'Fontsize',14);   
    end
    %***********
    subplot('position',[0.625 0.10 0.325 0.80]);
    smeanrf = imgaussfilt(meanrf,0.5);
    imagesc(Zx,Zy,smeanrf,[mino maxo]); hold on;
    plot([-15,15],[0,0],'k:','Linewidth',2);
    plot([0,0],[-15,15],'k:','Linewidth',2);
    plot([0,14.5],[14,14],'k-','Linewidth',2,'Color',[0.99,0.99,0.99]);
    ht = text(3,11,'15 dva','Fontsize',18,'Color',[0.99,0.99,0.99]);
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
    hc2 = text(Nx*1.7,Ny*0.6,'Firing Rate (sp/s)','Fontsize',16,'Rotation',90);
    %************

return;

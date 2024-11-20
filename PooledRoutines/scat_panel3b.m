function scat_panel3(xx,ndepth,nduration,narrow,broad,ncolo,bcolo,msize,vx,correlated,topname,lines,lamsmo,xtick,xtext,linep)
        
    xmin = min(vx);
    xmax = max(vx);
    
    grps = cell(2,3);
    lay = cell(1,1);
    cla = cell(1,1);
    obs = [];
    obcnt = 1;
    
    for layk = 1:3    
        
        subplot('Position',[xx (0.3+((layk-1)*0.22)) 0.24 0.18]);
        %******
        if (layk == 1)
           zz = find( ndepth(narrow) < 0);
           inarrow = narrow(zz);
           zz = find( ndepth(broad) < 0);
           ibroad = broad(zz);
           layname = '   Deep   ';
        end
        if (layk == 2)
           zz = find( (ndepth(narrow) < 300) & (ndepth(narrow) > 0) );
           inarrow = narrow(zz);
           zz = find( (ndepth(broad) > 300) & (ndepth(broad) > 0) );
           ibroad = broad(zz);
           layname = '   Input   ';
        end
        if (layk == 3)
           zz = find( ndepth(narrow) > 300);
           inarrow = narrow(zz);
           zz = find( ndepth(broad) > 300);
           ibroad = broad(zz);
           layname = 'Superficial';
        end
        grps{1,layk} = nduration(inarrow);
        grps{2,layk} = nduration(ibroad);
        %*********** code data into anovan format
        for ii = 1:length(inarrow)
            lay{1,obcnt} = sprintf('lay%d',layk);
            cla{1,obcnt} = 'cla1';
            obs(1,obcnt) = nduration(inarrow(ii));
            obcnt = obcnt + 1;
        end
        for ii = 1:length(ibroad)
            lay{1,obcnt} = sprintf('lay%d',layk);
            cla{1,obcnt} = 'cla2';
            obs(1,obcnt) = nduration(ibroad(ii));
            obcnt = obcnt + 1;
        end
        %*******  
        
        mhist = 0.50;
        minp = -0.005;
        aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
        bb = [minp,mhist,mhist,minp,minp];
        fill(aa,bb,[1,1,1],'FaceAlpha',1,'Linestyle','none'); hold on;
        plot([min(vx),max(vx)],[minp,minp],'k-','LineWidth',1);
        %plot([min(vx),max(vx)],[mhist,mhist],'k-','LineWidth',1);
        plot([min(vx),min(vx)],[minp,mhist],'k-','LineWidth',1);
        plot([max(vx),max(vx)],[minp,mhist],'k-','LineWidth',1);
        %*******
     
        nyx = hist(nduration(inarrow),vx);
        nyx = nyx / sum(nyx);
        byx = hist(nduration(ibroad),vx);
        byx = byx / sum(byx);
        bar(vx,nyx,1,'Linestyle','none','FaceColor',ncolo,'FaceAlpha',0.4); hold on;
        bar(vx,byx,1,'Linestyle','none','FaceColor',bcolo,'FaceAlpha',0.2);
        axis tight;
        V = axis;
        axis([xmin xmax V(3) (1.1*V(4))]);
        V = axis;
        if (correlated ~= 2) % && (correlated ~= 0)
            if (~correlated)
                xmen = 10.0;  % narrow vs broad thresh
            end
           % plot([xmen,xmen],[V(3),(1.1*V(4))],'k-');
        end
        %***
        if ~isempty(narrow)
           med1 = nanmedian(nduration(inarrow)); 
           if (correlated == 2)
            mod1 = med1;
           else
            mod1 = 100*(((1+med1)/(1-med1))-1);
           end
           p1 = signrank(nduration(inarrow));
        end
        %***
        if ~isempty(broad)
           med2 = nanmedian(nduration(ibroad)); 
           if (correlated == 2)
            mod2 = med2;
           else
            mod2 = 100*(((1+med2)/(1-med2))-1);
           end
           p2 = signrank(nduration(ibroad));
        end
        %*****
        if ~isempty(inarrow) && ~isempty(ibroad)
          p = ranksum(nduration(inarrow),nduration(ibroad));
          pa = signrank(nduration);
        else
          p = NaN;
          pa = NaN;
        end
        axis off;
        
        if (1) %(xx < 0.2)
          text((min(vx)-0.10*(max(vx)-min(vx))),0.0,sprintf('%3.1f',0.0),'Fontsize',12);
          text((min(vx)-0.10*(max(vx)-min(vx))),0.2,sprintf('%3.1f',0.2),'Fontsize',12);
          text((min(vx)-0.10*(max(vx)-min(vx))),0.4,sprintf('%3.1f',0.4),'Fontsize',12);
          text((min(vx)-0.17*(max(vx)-min(vx))),0.025,'Probability','Fontsize',14,'Rotation',90);
        end
        plot([min(vx),(min(vx)+0.03*(max(vx)-min(vx)))],[0.2,0.2],'k-','Linewidth',1);
        plot([min(vx),(min(vx)+0.03*(max(vx)-min(vx)))],[0.4,0.4],'k-','Linewidth',1);
  
        if (1)
          dm = max(vx)-min(vx);  
          xm = 0.05*dm;
          if (1)
             plot([0,0],[0,mhist],'k-','Color',[0.0,0.0,0.0]);
          end 
          if (layk == 1)
            for tk = 1:length(xtick)
             text(xtick(tk)-xm,-0.06,xtext{tk},'Fontsize',12);
             plot([xtick(1),xtick(1)],[-0.01,-0.02],'k-');
            end
            if (1)
             text(min(vx)+(dm/8),-0.15,topname,'Fontsize',14);
             plot([0,0],[0,mhist],'k-','Color',[0.0,0.0,0.0]);
            else 
             text(min(vx)+(dm/4),-0.15,topname,'Fontsize',14);
            end
          end
          text(min(vx)+(1.5*dm/4),1.1*mhist,layname,'Fontsize',14);
       
          %********
          if ~isempty(inarrow)
             tx = text(min(vx)+(3.3*dm/5),0.75*mhist,sprintf(' Median %3.1f',mod1),'Fontsize',12,'Color',0.7*ncolo);          
             if (p < 0.05)
                 set(tx,'FontWeight','Bold');
             end
          end
          if ~isempty(ibroad)
             tx = text(min(vx)+(3.3*dm/5),0.90*mhist,sprintf(' Median %3.1f',mod2),'Fontsize',12,'Color',0.7*bcolo);          
             if (p < 0.05)
                 set(tx,'FontWeight','Bold');
             end
          end
          if ~isempty(inarrow) & ~isempty(ibroad)
             tx = text(min(vx)+(3.4*dm/5),0.60*mhist,sprintf(' P diff %5.3f',p),'Fontsize',12,'Color',[0,0,0]);                          
             if (p < 0.05)
                 set(tx,'FontWeight','Bold');
             end
          end
        end   
    end
    %****************
    [~,~,stats] = anovan(obs,{lay cla},'model','interaction','varnames',{'lay','cla'});
    % input('check');
    %****************
           
return;



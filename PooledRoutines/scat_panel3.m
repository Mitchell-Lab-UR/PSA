function scat_panel3(xx,ndepth,nduration,narrow,broad,ncolo,bcolo,msize,vx,correlated,topname,lines,lamsmo,xtick,xtext,linep)
        %************
        subplot('Position',[xx 0.30 0.25 0.625]); hold off;
        %***** lay down a white box
        mndepth = min(ndepth)-50;
        mxdepth = max(ndepth)+50;
        aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
        bb = [mndepth,mxdepth,mxdepth,mndepth,mndepth];
        fill(aa,bb,[1,1,1],'FaceAlpha',1,'Linestyle','none'); hold on;
        if (linep==3)
         if (length(lines)>2)  % gray out intermediate depth zones
          aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
          bb = [lines(1),lines(3),lines(3),lines(1),lines(1)];
          fill(aa,bb,[0.92,0.92,0.92],'FaceAlpha',1,'Linestyle','none'); hold on;
          aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
          bb = [lines(4),lines(2),lines(2),lines(4),lines(4)];
          fill(aa,bb,[0.92,0.92,0.92],'FaceAlpha',1,'Linestyle','none'); hold on;    
         end
        end
        plot([min(vx),max(vx)],[mndepth,mndepth],'k-','LineWidth',1);
        plot([min(vx),max(vx)],[mxdepth,mxdepth],'k-','LineWidth',1);
        plot([min(vx),min(vx)],[mndepth,mxdepth],'k-','LineWidth',1);
        plot([max(vx),max(vx)],[mndepth,mxdepth],'k-','LineWidth',1);
        %******
        if (correlated)
          % plot(nduration,ndepth,'k.','Markersize',(msize/2)); hold on;
          %******
          xd = 0.0125;
          yd = 5;
          ang = 0:0.1:(2*pi);
          xdd = xd * cos(ang);
          ydd = yd * sin(ang);
          %********
          if ~isempty(narrow)
             for di = 1:length(narrow)
                xi = nduration(narrow(di));
                yi = ndepth(narrow(di));
                aa = [(xi-xd),(xi+xd),(xi+xd),(xi-xd),(xi-xd)];
                bb = [(yi-yd),(yi-yd),(yi+yd),(yi+yd),(yi-yd)];
                if (0) % linep == 2)
                   plot(aa,bb,'k-','Color',ncolo); hold on;        
                else  
                   for k = 1:10
                      sp = k;
                      fill((sp*xdd+xi),(sp*ydd+yi),ncolo,'FaceAlpha',0.05,'Linestyle','none'); hold on;
                   end
                end
             end
             for di = 1:length(narrow)
                xi = nduration(narrow(di));
                yi = ndepth(narrow(di));
                aa = [(xi-xd),(xi+xd),(xi+xd),(xi-xd),(xi-xd)];
                bb = [(yi-yd),(yi-yd),(yi+yd),(yi+yd),(yi-yd)];
                if (0) % (linep == 2)
                   plot(aa,bb,'k-','Color',ncolo); hold on;        
                else  
                   sp = 1;
                   fill((sp*xdd+xi),(sp*ydd+yi),ncolo*0.5,'FaceAlpha',1.0,'Linestyle','none'); hold on;
                end
             end   
             %****** add regression line
             if (linep == 1)
               zx = nduration(narrow);
               zy = ndepth(narrow);
               b = gmregress(zy,zx);  % Model II regression line
               zyy = mndepth:20:mxdepth;
               zxx = b(2)*zyy + b(1);
               plot(zxx,zyy,'k-','Color',ncolo*1.0,'LineWidth',2);
             end
             %*********
          end
          if ~isempty(broad) 
             for di = 1:length(broad)
                xi = nduration(broad(di));
                yi = ndepth(broad(di));
                aa = [(xi-xd),(xi+xd),(xi+xd),(xi-xd),(xi-xd)];
                bb = [(yi-yd),(yi-yd),(yi+yd),(yi+yd),(yi-yd)];
                if (0) % linep == 2)
                   plot(aa,bb,'k-','Color',bcolo); hold on;        
                else  
                   for k = 1:10
                      sp = k;
                      fill((sp*xdd+xi),(sp*ydd+yi),bcolo,'FaceAlpha',0.05,'Linestyle','none'); hold on;
                   end
                end
             end
             for di = 1:length(broad)
                xi = nduration(broad(di));
                yi = ndepth(broad(di));
                aa = [(xi-xd),(xi+xd),(xi+xd),(xi-xd),(xi-xd)];
                bb = [(yi-yd),(yi-yd),(yi+yd),(yi+yd),(yi-yd)];
                if (0) % linep == 2)
                   plot(aa,bb,'k-','Color',ncolo); hold on;        
                else  
                   sp = 1;
                   fill((sp*xdd+xi),(sp*ydd+yi),bcolo*0.5,'FaceAlpha',1.0,'Linestyle','none'); hold on;
                end
             end
             %****** add regression line
             if (linep == 1)
               zx = nduration(broad);
               zy = ndepth(broad);
               % XX = [zy ones(size(zy))];
               % YY = zx;
               % b = regress(YY,XX);
               b = gmregress(zy,zx);
               zyy = mndepth:20:mxdepth;
               % zxx = b(1)*zyy + b(2);
               zxx = b(2)*zyy + b(1);
               plot(zxx,zyy,'k-','Color',bcolo*1.0,'LineWidth',2);
             end
          end   
          %plot(nduration(narrow),ndepth(narrow),'ko','Markersize',msize/2,'Color',ncolo); hold on;
          %plot(nduration(broad),ndepth(broad),'ko','Markersize',msize/2,'Color',bcolo); hold on;
        else  
          plot(nduration,ndepth,'ko','Markersize',msize); hold on;  % show in black only is cluster plot
         % plot(nduration(narrow),ndepth(narrow),'ko','Markersize',msize,'Color',ncolo); hold on;
         % plot(nduration(broad),ndepth(broad),'ko','Markersize',msize,'Color',bcolo); hold on;     
        end
        %******* compute rolling mean over depth
        if (correlated)
            if (linep == 2)  % sliding depth window
              dmin = mndepth;
              dmax = mxdepth;
              dbin = lamsmo;
              minum = 2;
              nuu = [];
              nsu = [];
              buu = [];
              bsu = [];
              auu = [];
              asu = [];
              xuu = [];
              % for dx = dmin:((dmax-dmin)/100):dmax
              for dx = dmin:(dbin/2):dmax
                dxa = dx-dbin;
                dxb = dx+dbin;
                xuu = [xuu ; dx];
                %******
                zz = find( (ndepth(narrow) > dxa) & (ndepth(narrow) < dxb));
                if length(zz) < minum
                    nuu = [nuu ; NaN];
                    nsu = [nsu ; NaN];
                else
                    nuu = [nuu ; nanmean(nduration(narrow(zz)))];
                    nsu = [nsu ; (nanstd(nduration(narrow(zz)))/sqrt(length(zz)))];
                end
                %*******
                zz = find( (ndepth(broad) > dxa) & (ndepth(broad) < dxb));
                if length(zz) < minum
                    buu = [buu ; NaN];
                    bsu = [bsu ; NaN];
                else
                    buu = [buu ; nanmean(nduration(broad(zz)))];
                    bsu = [bsu ; (nanstd(nduration(broad(zz)))/sqrt(length(zz)))];
                end
                %*******
                zz = find( (ndepth > dxa) & (ndepth < dxb));
                if length(zz) < minum
                    auu = [auu ; NaN];
                    asu = [asu ; NaN];
                else
                    auu = [auu ; nanmean(nduration(zz))];
                    asu = [asu ; (nanstd(nduration(zz))/sqrt(length(zz)))];
                end
              end
              %******** 
              iz = find( ~isnan(nuu));
              aa = [xuu(iz) ; flipud(xuu(iz))];
              bb = [(nuu(iz)+(2*nsu(iz))) ; flipud( nuu(iz)-(2*nsu(iz)))];
              fill(bb,aa,ncolo*0.0,'FaceAlpha',0.25,'Linestyle','none');
              plot(nuu(iz),xuu(iz),'k-','Color',ncolo*0.0,'Linewidth',3);
              %plot(nuu+(2*nsu),xuu,'k-','Color',ncolo,'Linewidth',1);
              %plot(nuu-(2*nsu),xuu,'k-','Color',ncolo,'Linewidth',1);
              %******
              iz = find( ~isnan(buu));
              aa = [xuu(iz) ; flipud(xuu(iz))];
              bb = [(buu(iz)+(2*bsu(iz))) ; flipud( buu(iz)-(2*bsu(iz)))];
              fill(bb,aa,bcolo*0.0,'FaceAlpha',0.25,'Linestyle','none');
              plot(buu,xuu,'k-','Color',bcolo*0.0,'Linewidth',3);
              % plot(buu+(2*bsu),xuu,'k-','Color',bcolo,'Linewidth',1);
              % plot(buu-(2*bsu),xuu,'k-','Color',bcolo,'Linewidth',1);
            end
            %********
            %******
            if (0)
              iz = find( ~isnan(auu));
              aa = [xuu(iz) ; flipud(xuu(iz))];
              bb = [(auu(iz)+(2*asu(iz))) ; flipud( auu(iz)-(2*asu(iz)))];
              fill(bb,aa,[0,0,0],'FaceAlpha',0.2,'Linestyle','none');
              plot(auu,xuu,'k-','Linewidth',2);
            end
            %********
        end
        
        axis tight;
        V = axis;
        xmin = min(vx);
        xmax = max(vx);
        axis([xmin xmax V(3) V(4)]);
        V = axis;
        xmen = 0.5*(xmin+xmax);
        if (linep ~= 0)
          if (correlated ~= 2) && (correlated ~= 0) 
             plot([xmen,xmen],[V(3),V(4)],'k-');
          end
        else
          if (correlated ~= 2) && (correlated ~= 0) 
             plot([xmen,xmen],[V(3),V(4)],'k-','Color',[0.2,0.2,0.2]);
          end       
        end
        if ~isempty(lines) % && (length(lines) <= 2)
          plot([V(1),V(2)],[lines(1),lines(1)],'k--'); %LineWidth,2);
          %**** hard coded for now
          if (1)
            plot([V(1),V(2)],[lines(2),lines(2)],'r--');     
          else
            plot([V(1),V(2)],[lines(2),lines(2)],'k:');
          end
        end
        %*********
        cotype = 'Spearman';   %middle
        % cotype = 'Kendall';  % most conservative
        % cotype = 'Pearson';  % not consertative
        
        if (correlated)
          if ~isempty(narrow)  
            bz = find( ~isnan(nduration(narrow)) & ~isnan(ndepth(narrow)) );
            [r1,p1] = corr(nduration(narrow(bz)),ndepth(narrow(bz)),'type',cotype);
          else
            r1 = 0; p1 = 0;  
          end
          if ~isempty(broad)
             bz = find( ~isnan(nduration(broad)) & ~isnan(ndepth(broad)) );
            [r2,p2] = corr(nduration(broad(bz)),ndepth(broad(bz)),'type',cotype);
          else
             r2 = 0; p2 = 1;
          end
          total = [narrow ; broad];
          tz = find( ~isnan(nduration(total)) & ~isnan(ndepth(total)) );
          [rt,pt] = corr(nduration(total(tz)),ndepth(total(tz)),'type',cotype);
          if (linep > -1)
            if ~isempty(narrow)
             text(V(1)+0.2*(V(2)-V(1)),V(4)+0.04*(V(4)-V(3)),sprintf('R %4.3f   p=%5.3f',r1,p1),'Color',ncolo,'Fontsize',18);
             text(0.35,750,sprintf('N=%d',length(narrow)),'Color',ncolo,'Fontsize',16);
            else
             text(V(1)+0.2*(V(2)-V(1)),V(4)+0.04*(V(4)-V(3)),sprintf('R %4.3f   p=%5.3f',r2,p2),'Color',bcolo,'Fontsize',18);
             text(-0.85,750,sprintf('N=%d',length(broad)),'Color',bcolo,'Fontsize',16);            
            end
          end
          % title(sprintf('N(%d)%4.2f(%5.3f);B(%d)%4.2f(%5.3f)',length(narrow),r1,p1,length(broad),r2,p2));
          %******
          %disp('   ');
          %disp(sprintf('%s : Depth correlations',topname));
          %disp(sprintf('Nw: %6.4f(%8.6f); Bd:%6.4f(%8.6f); Tot:%6.4f(%8.6f)',r1,p1,r2,p2,rt,pt));
         
          %******** plot laminar Y-axis
          if (xx < 0.2)
            plot([min(vx),min(vx)+0.04*(max(vx)-min(vx))],[lines(1),lines(1)],'k-');
            plot([min(vx),min(vx)+0.04*(max(vx)-min(vx))],[lines(2),lines(2)],'k-');
            text((min(vx)-0.16*(max(vx)-min(vx))),lines(1),sprintf('%3d',floor(lines(1))),'Fontsize',12);
            text((min(vx)-0.16*(max(vx)-min(vx))),lines(2),sprintf('%3d',floor(lines(2))),'Fontsize',12);
          end
          
          %******
          % axis off;
        else
          %**** put in vertical tick marks
          if (xx < 0.2)
            plot([min(vx),min(vx)+0.04*(max(vx)-min(vx))],[lines(1),lines(1)],'k-');
            plot([min(vx),min(vx)+0.04*(max(vx)-min(vx))],[lines(2),lines(2)],'k-');
            text((min(vx)-0.16*(max(vx)-min(vx))),lines(1),sprintf('%3d',floor(lines(1))),'Fontsize',12);
            text((min(vx)-0.16*(max(vx)-min(vx))),lines(2),sprintf('%3d',floor(lines(2))),'Fontsize',12);
            if (length(lines)>2)
             plot([min(vx),min(vx)+0.04*(max(vx)-min(vx))],[lines(3),lines(3)],'k-');
             plot([min(vx),min(vx)+0.04*(max(vx)-min(vx))],[lines(4),lines(4)],'k-');
             text((min(vx)-0.16*(max(vx)-min(vx))),lines(3),sprintf('%3d',floor(lines(3))),'Fontsize',12);
             text((min(vx)-0.16*(max(vx)-min(vx))),lines(4),sprintf('%3d',floor(lines(4))),'Fontsize',12);
            end
            text((min(vx)-0.30*(max(vx)-min(vx))),lines(2)-50,'Cortical Depth (um)','Fontsize',18,'Rotation',90);
            title(sprintf('Narrow(N=%d) Broad(N=%d)',length(narrow),length(broad)));
          end
        end
        axis off;
        if (xx < 0.2)
          text((min(vx)-0.30*(max(vx)-min(vx))),lines(2)-50,'Cortical Depth (um)','Fontsize',18,'Rotation',90);
        end
        
        subplot('Position',[xx 0.1 0.25 0.20]);
        mhist = 0.35; %35;
        minp = -0.005;
        aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
        bb = [minp,mhist,mhist,minp,minp];
        fill(aa,bb,[1,1,1],'FaceAlpha',1,'Linestyle','none'); hold on;
        plot([min(vx),max(vx)],[minp,minp],'k-','LineWidth',1);
        %plot([min(vx),max(vx)],[mhist,mhist],'k-','LineWidth',1);
        plot([min(vx),min(vx)],[minp,mhist],'k-','LineWidth',1);
        plot([max(vx),max(vx)],[minp,mhist],'k-','LineWidth',1);
        %*******
        if ~correlated
           yx = hist(nduration,vx);
           nyx = hist(nduration(narrow),vx);
           byx = hist(nduration(broad),vx);
           %******
           if (0)
             omax = max([max(nyx),max(byx),max(yx)]);
             yx = (yx * (0.8*mhist/omax));
             nyx = (nyx * (0.8*mhist/omax));
             byx = (byx * (0.8*mhist/omax));
           else
             omax = sum(yx);
             yx = 2 * yx / omax;
             nyx = 2 * nyx / omax;
             byx = 2 * byx / omax;
           end
           %******
           bar(vx,yx,1,'Linestyle','none','FaceColor',[0,0,0],'FaceAlpha',0.3); hold on;
           %****** superimpose colors on top
           bva = 0.6;
           bar(vx,nyx,1,'Linestyle','none','FaceColor',...
                      [bva,bva,bva]+(1-bva)*ncolo,'FaceAlpha',1); hold on;
           bar(vx,byx,1,'Linestyle','none','FaceColor',...
                      [bva,bva,bva]+(1-bva)*bcolo,'FaceAlpha',1);
        else
          nyx = hist(nduration(narrow),vx);
          nyx = nyx / sum(nyx);
          byx = hist(nduration(broad),vx);
          byx = byx / sum(byx);
          % plot(vx,nyx,'k.-','Color',ncolo,'Linewidth',2); hold on;
          % plot(vx,byx,'k.-','Color',bcolo,'Linewidth',2);
          bar(vx,nyx,1,'Linestyle','none','FaceColor',ncolo,'FaceAlpha',0.4); hold on;
          bar(vx,byx,1,'Linestyle','none','FaceColor',bcolo,'FaceAlpha',0.2);
          % plot([xmen,xmen],[V(3),V(4)],'k--');
        end
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
           med1 = nanmedian(nduration(narrow)); 
           if (correlated == 2)
            mod1 = med1;
           else
            mod1 = 100*(((1+med1)/(1-med1))-1);
           end
           p1 = signrank(nduration(narrow));
        end
        %***
        if ~isempty(broad)
           med2 = nanmedian(nduration(broad)); 
           if (correlated == 2)
            mod2 = med2;
           else
            mod2 = 100*(((1+med2)/(1-med2))-1);
           end
           p2 = signrank(nduration(broad));
        end
        %*****
        if ~isempty(narrow) && ~isempty(broad)
          p = ranksum(nduration(narrow),nduration(broad));
          pa = signrank(nduration);
        else
          p = NaN;
          pa = NaN;
        end
        if (correlated)
            % title(sprintf('N:%5.2f,B:%5.2f(p=%4.3f)',... %A:%5.2f(p=%4.3f)',...
            %       mod1,mod2,p)); %moda,pa));
            % text(xmin+0.65*(xmax-xmin),0.90*V(4),sprintf('N:%5.2f(p=%4.3f)',mod1,p1),'Fontsize',8,'Color',ncolo);
            % text(xmin+0.65*(xmax-xmin),0.80*V(4),sprintf('B:%5.2f(p=%4.3f)',mod2,p2),'Fontsize',8,'Color',bcolo);   
        
            %******* display stats on the distributions narrow vs broad
            disp(sprintf('%s : Narrow vs Broad Without Layer',topname));
            if ~isempty(narrow)
              disp(sprintf('Narrow: Mod %5.3f (p=%8.20f)',mod1,p1));
              disp(p1);
            end
            if ~isempty(broad)
              disp(sprintf('Broad : Mod %5.3f (p=%8.20f)',mod2,p2));
              disp(p2);
            end
            if ~isempty(narrow) & ~isempty(broad)
              disp(sprintf('Difference Comparison (p=%8.20f)',p));
              disp(p);
              disp(sprintf('Total:  Mod %5.3f (p=%8.2pf)',moda,pa));
              disp(pa);
            end
            disp('   ');
            %*********
        end
        axis off;
        
        if (xx < 0.2)
          text((min(vx)-0.12*(max(vx)-min(vx))),0.0,sprintf('%3.1f',0.0),'Fontsize',12);
          text((min(vx)-0.12*(max(vx)-min(vx))),0.1,sprintf('%3.1f',0.1),'Fontsize',12);
          text((min(vx)-0.12*(max(vx)-min(vx))),0.2,sprintf('%3.1f',0.2),'Fontsize',12);
          text((min(vx)-0.25*(max(vx)-min(vx))),0.025,'Probability','Fontsize',18,'Rotation',90);
        end
        plot([min(vx),(min(vx)+0.03*(max(vx)-min(vx)))],[0.1,0.1],'k-','Linewidth',1);
        plot([min(vx),(min(vx)+0.03*(max(vx)-min(vx)))],[0.2,0.2],'k-','Linewidth',1);
         
        if ~correlated
          text((min(vx)-0.12*(max(vx)-min(vx))),0.0,sprintf('%3.1f',0.0),'Fontsize',12);
          text((min(vx)-0.12*(max(vx)-min(vx))),0.2,sprintf('%3.1f',0.2),'Fontsize',12);
          text((min(vx)-0.20*(max(vx)-min(vx))),0.035,'Probability','Fontsize',18,'Rotation',90);
          dm = max(vx)-min(vx);
          xm = 0.05*dm;
          for tk = 1:length(xtick)
             text(xtick(tk)-xm,-0.04,xtext{tk},'Fontsize',12);
             plot([xtick(1),xtick(1)],[-0.01,-0.02],'k-');
          end
          text(0.30,-0.10,topname,'Fontsize',14);
          %***** show dip test
          text(min(vx)+(1*dm/12),0.95*mhist,sprintf('Narrow: N=%d',length(narrow)),'Fontsize',10,'Color',ncolo);          
          text(min(vx)+(1*dm/12),0.85*mhist,sprintf(' Broad: N=%d',length(broad)),'Fontsize',10,'Color',bcolo);     
          hp = HartigansDipTest(nduration);
          text(min(vx)+(1*dm/12),0.75*mhist,sprintf(' Dip-test(p=%4.3f)',hp),'Fontsize',10); 
        else
          dm = max(vx)-min(vx);  
          xm = 0.05*dm;
          for tk = 1:length(xtick)
             text(xtick(tk)-xm,-0.04,xtext{tk},'Fontsize',12);
             plot([xtick(1),xtick(1)],[-0.01,-0.02],'k-');
          end
          if (length(topname)>20)
             text(min(vx)+(dm/16),-0.10,topname,'Fontsize',14);
             plot([0,0],[0,mhist],'k-','Color',[0.2,0.2,0.2]);
          else 
             text(min(vx)+(dm/4),-0.10,topname,'Fontsize',14);
          end
          %********
          if ~isempty(narrow)
             text(min(vx)+(3*dm/5),0.90*mhist,sprintf('Median %3.1f',mod1),'Fontsize',12,'Color',0.7*ncolo);          
             text(min(vx)+(3*dm/5),0.75*mhist,sprintf('  p = %6.4f',p1),'Fontsize',12,'Color',0.7*ncolo);          
          end
          if ~isempty(broad)
             text(min(vx)+(3*dm/5),0.90*mhist,sprintf('Median %3.1f',mod2),'Fontsize',12,'Color',0.7*bcolo);          
             text(min(vx)+(3*dm/5),0.75*mhist,sprintf('  p = %6.4f',p2),'Fontsize',12,'Color',0.7*bcolo);                      
          end
         end
        
        %******************************
        %******** Break out all stats by layers
        for lay = 1:3
           %********* 
           if (lay == 1)
               inarrow = narrow( find( (ndepth(narrow) > 300) ));
               ibroad = broad( find( (ndepth(broad) > 300) ));
           end
           if (lay == 2)
               inarrow = narrow( find( (ndepth(narrow) <= 300) & ...
                                       (ndepth(narrow) >= 0)) );
               ibroad = broad( find( (ndepth(broad) <= 300) & ...
                                     (ndepth(broad) >= 0)) );
           end
           if (lay == 3)
               inarrow = narrow( find( (ndepth(narrow) < 0) ));
               ibroad = broad( find( (ndepth(broad) < 0) ));
           end
           itotal = [inarrow ; ibroad];
           %*******
           med1 = nanmedian(nduration(inarrow)); 
           med2 = nanmedian(nduration(ibroad));
           meda = nanmedian(nduration(itotal));
           if (correlated == 2)
             mod1 = med1;
             mod2 = med2;
             moda = meda;
           else
             mod1 = 100*(((1+med1)/(1-med1))-1);
             mod2 = 100*(((1+med2)/(1-med2))-1);
             moda = 100*(((1+meda)/(1-meda))-1);
           end
           if (0)
               p = ranksum(nduration(inarrow),nduration(ibroad));
               p1 = signrank(nduration(inarrow));
               p2 = signrank(nduration(ibroad));
               pa = signrank(nduration(itotal));
               %********************* 
               disp(sprintf('>>>> RESULTS BY LAYER %d, %s',lay,topname));
               disp(sprintf('>>>> N: Narrow(%d) Broad(%d) Tot(%d)',...
                         length(inarrow),length(ibroad),length(itotal)));
               disp(sprintf('>>>> Narrow: Mod %5.3f (p=%8.20f)',mod1,p1));
               disp(p1);
               disp(sprintf('>>>> Broad : Mod %5.3f (p=%8.20f)',mod2,p2));
               disp(p2);
               disp(sprintf('>>>> Difference Comparison (p=%8.20f)',p));
               disp(p);
               disp(sprintf('>>>> Total:  Mod %5.3f (p=%8.20f)',moda,pa));
               disp(pa);
               disp('   ');
           end
           %******************** 
        end
        %****************
         
           
return;



function scat_panel2(xx,ndepth,nduration,narrow,broad,ncolo,bcolo,msize,vx,correlated,topname,lines,lamsmo,xtick,xtext)
        %************
        subplot('Position',[xx 0.30 0.20 0.65]); hold off;
        %***** lay down a white box
        aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
        bb = [min(ndepth),max(ndepth),max(ndepth),min(ndepth),min(ndepth)];
        fill(aa,bb,[1,1,1],'FaceAlpha',1,'Linestyle','none'); hold on;
        if (length(lines)>2)  % gray out intermediate depth zones
          aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
          bb = [lines(1),lines(3),lines(3),lines(1),lines(1)];
          fill(aa,bb,[0.92,0.92,0.92],'FaceAlpha',1,'Linestyle','none'); hold on;
          aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
          bb = [lines(4),lines(2),lines(2),lines(4),lines(4)];
          fill(aa,bb,[0.92,0.92,0.92],'FaceAlpha',1,'Linestyle','none'); hold on;    
        end
        plot([min(vx),max(vx)],[min(ndepth),min(ndepth)],'k-','LineWidth',1);
        plot([min(vx),max(vx)],[max(ndepth),max(ndepth)],'k-','LineWidth',1);
        plot([min(vx),min(vx)],[min(ndepth),max(ndepth)],'k-','LineWidth',1);
        plot([max(vx),max(vx)],[min(ndepth),max(ndepth)],'k-','LineWidth',1);
        %******
        if (correlated)
          % plot(nduration,ndepth,'k.','Markersize',(msize/2)); hold on;
          plot(nduration(narrow),ndepth(narrow),'ko','Markersize',msize/2,'Color',ncolo); hold on;
          plot(nduration(broad),ndepth(broad),'k.','Markersize',msize,'Color',bcolo); hold on;
        else  
          plot(nduration,ndepth,'ko','Markersize',msize); hold on;  % show in black only is cluster plot
         % plot(nduration(narrow),ndepth(narrow),'ko','Markersize',msize,'Color',ncolo); hold on;
         % plot(nduration(broad),ndepth(broad),'ko','Markersize',msize,'Color',bcolo); hold on;     
        end
        %******* compute rolling mean over depth
        if (correlated)
            dmin = min(ndepth);
            dmax = max(ndepth);
            dbin = lamsmo;
            minum = 2;
            nuu = [];
            nsu = [];
            buu = [];
            bsu = [];
            auu = [];
            asu = [];
            xuu = [];
            for dx = dmin:((dmax-dmin)/100):dmax
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
            fill(bb,aa,ncolo,'FaceAlpha',0.2,'Linestyle','none');
            plot(nuu(iz),xuu(iz),'k-','Color',ncolo,'Linewidth',2);
            %plot(nuu+(2*nsu),xuu,'k-','Color',ncolo,'Linewidth',1);
            %plot(nuu-(2*nsu),xuu,'k-','Color',ncolo,'Linewidth',1);
            %******
            iz = find( ~isnan(buu));
            aa = [xuu(iz) ; flipud(xuu(iz))];
            bb = [(buu(iz)+(2*bsu(iz))) ; flipud( buu(iz)-(2*bsu(iz)))];
            fill(bb,aa,bcolo,'FaceAlpha',0.2,'Linestyle','none');
            plot(buu,xuu,'k-','Color',bcolo,'Linewidth',2);
            % plot(buu+(2*bsu),xuu,'k-','Color',bcolo,'Linewidth',1);
            % plot(buu-(2*bsu),xuu,'k-','Color',bcolo,'Linewidth',1);
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
        if (correlated ~= 2) && (correlated ~= 0)
           plot([xmen,xmen],[V(3),V(4)],'k-');
        end
        if ~isempty(lines) && (length(lines) <= 2)
          plot([V(1),V(2)],[lines(1),lines(1)],'r-'); %LineWidth,2);
          %**** hard coded for now
          if (1)
            plot([V(1),V(2)],[lines(2),lines(2)],'k-');     
          else
            plot([V(1),V(2)],[lines(2),lines(2)],'k:');
          end
        end
        %*********
        if (correlated)
          if ~isempty(narrow)  
            bz = find( ~isnan(nduration(narrow)) & ~isnan(ndepth(narrow)) );
            [r1,p1] = corr(nduration(narrow(bz)),ndepth(narrow(bz)),'type','Spearman');
          else
            r1 = 0; p1 = 0;  
          end
          if ~isempty(broad)
             bz = find( ~isnan(nduration(broad)) & ~isnan(ndepth(broad)) );
            [r2,p2] = corr(nduration(broad(bz)),ndepth(broad(bz)),'type','Spearman');
          else
             r2 = 0; p2 = 1;
          end
          total = [narrow ; broad];
          tz = find( ~isnan(nduration(total)) & ~isnan(ndepth(total)) );
          [rt,pt] = corr(nduration(total(tz)),ndepth(total(tz)),'type','Spearman');
          title(sprintf('Nw:%4.2f(%5.3f);Bd:%4.2f(%5.3f)',r1,p1,r2,p2));
          %******
          disp('   ');
          disp(sprintf('%s : Depth correlations',topname));
          disp(sprintf('Nw: %6.4f(%8.6f); Bd:%6.4f(%8.6f); Tot:%6.4f(%8.6f)',r1,p1,r2,p2,rt,pt));
          %******
          % axis off;
        else
          %**** put in vertical tick marks
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
          text((min(vx)-0.30*(max(vx)-min(vx))),lines(2)-50,'Cortical Depth (um)','Fontsize',16,'Rotation',90);
          title(sprintf('Narrow(N=%d) Broad(N=%d)',length(narrow),length(broad)));
        end
        axis off;
        
        
        subplot('Position',[xx 0.1 0.20 0.20]);
        mhist = 0.35;
        minp = -0.005;
        aa = [min(vx),min(vx),max(vx),max(vx),min(vx)];
        bb = [minp,mhist,mhist,minp,minp];
        fill(aa,bb,[1,1,1],'FaceAlpha',1,'Linestyle','none'); hold on;
        plot([min(vx),max(vx)],[minp,minp],'k-','LineWidth',1);
        plot([min(vx),max(vx)],[mhist,mhist],'k-','LineWidth',1);
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
        med1 = nanmedian(nduration(narrow)); 
        med2 = nanmedian(nduration(broad));
        meda = nanmedian(nduration);
        if (correlated == 2)
            mod1 = med1;
            mod2 = med2;
            moda = meda;
        else
           mod1 = 100*(((1+med1)/(1-med1))-1);
           mod2 = 100*(((1+med2)/(1-med2))-1);
           moda = 100*(((1+meda)/(1-meda))-1);
        end
        if ~isempty(narrow) && ~isempty(broad)
          p = ranksum(nduration(narrow),nduration(broad));
          p1 = signrank(nduration(narrow));
          p2 = signrank(nduration(broad));
          pa = signrank(nduration);
        else
          p = NaN;
          p1 = NaN;
          p2 = NaN;
          pa = NaN;
        end
        if (correlated)
            % title(sprintf('N:%5.2f,B:%5.2f(p=%4.3f)',... %A:%5.2f(p=%4.3f)',...
            %       mod1,mod2,p)); %moda,pa));
            % text(xmin+0.65*(xmax-xmin),0.90*V(4),sprintf('N:%5.2f(p=%4.3f)',mod1,p1),'Fontsize',8,'Color',ncolo);
            % text(xmin+0.65*(xmax-xmin),0.80*V(4),sprintf('B:%5.2f(p=%4.3f)',mod2,p2),'Fontsize',8,'Color',bcolo);   
        
            %******* display stats on the distributions narrow vs broad
            disp(sprintf('%s : Narrow vs Broad Without Layer',topname));
            disp(sprintf('Narrow: Mod %5.3f (p=%8.20f)',mod1,p1));
            disp(p1);
            disp(sprintf('Broad : Mod %5.3f (p=%8.20f)',mod2,p2));
            disp(p2);
            disp(sprintf('Difference Comparison (p=%8.20f)',p));
            disp(p);
            disp(sprintf('Total:  Mod %5.3f (p=%8.2pf)',moda,pa));
            disp(pa);
            disp('   ');
            %*********
        end
        axis off;
        if ~correlated
          text((min(vx)-0.12*(max(vx)-min(vx))),0.0,sprintf('%3.1f',0.0),'Fontsize',12);
          text((min(vx)-0.12*(max(vx)-min(vx))),0.2,sprintf('%3.1f',0.2),'Fontsize',12);
          text((min(vx)-0.25*(max(vx)-min(vx))),0.065,'Probability','Fontsize',16,'Rotation',90);
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
          else 
             text(min(vx)+(dm/4),-0.10,topname,'Fontsize',14);
          end
          %********
          %if (max(vx) < 20)
          %    mod1 = 10^mod1;
          %    mod2 = 10^mod2;
          %end
          text(min(vx)+(3*dm/5),0.95*mhist,sprintf('Nw Med: %4.2f',mod1),'Fontsize',10,'Color',ncolo);          
          text(min(vx)+(3*dm/5),0.85*mhist,sprintf(' Bd Med: %4.2f',mod2),'Fontsize',10,'Color',bcolo);          
          text(min(vx)+(3*dm/5),0.75*mhist,sprintf('   (p=%5.4f)',p),'Fontsize',10);
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



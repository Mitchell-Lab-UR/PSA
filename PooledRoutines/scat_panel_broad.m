function scat_panel_broad(xx,ndepth,nduration,narrow,broad,ncolo,bcolo,msize,vx,correlated,topname)
        %************
        subplot('Position',[xx 0.35 0.20 0.60]);
        if (correlated)
          %plot(nduration,ndepth,'k.','Markersize',(msize/2)); hold on;
          %plot(nduration(narrow),ndepth(narrow),'k.','Markersize',msize,'Color',ncolo); hold on;
          plot(nduration(broad),ndepth(broad),'k.','Markersize',msize,'Color',bcolo); hold on;
        else  
          %plot(nduration,ndepth,'ko','Markersize',msize); hold on;
          %plot(nduration(narrow),ndepth(narrow),'ko','Markersize',msize,'Color',ncolo); hold on;
          plot(nduration(broad),ndepth(broad),'ko','Markersize',msize,'Color',bcolo); hold on;     
        end
        %******* compute rolling mean over depth
        if (correlated)
            dmin = min(ndepth);
            dmax = max(ndepth);
            dbin = 100;
            %minum = 10;
            minum = 5;
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
%             fill(bb,aa,ncolo,'FaceAlpha',0.2,'Linestyle','none');
%             plot(nuu(iz),xuu(iz),'k-','Color',ncolo,'Linewidth',2);
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
        
%         axis tight;
%         V = axis;
%         xmin = min(vx);
%         xmax = max(vx);
%         axis([xmin xmax V(3) V(4)]);
%         V = axis;
        axis([-0.55 0.55 -400 800]); 
%         xmen = 0.5*(xmin+xmax);
        if (correlated ~= 2) && (correlated ~= 0)
%            plot([xmen,xmen],[V(3),V(4)],'k-');
             plot([0,0],[-400 800],'k-');
        end
%         plot([V(1),V(2)],[0,0],'r--');
          plot([-0.55 0.55],[0,0],'r--');
        %**** hard coded for now
%         plot([V(1),V(2)],[300,300],'k--');
         plot([-0.55 0.55],[300,300],'k--');
        %*********
        if (correlated)
          bz = find( ~isnan(nduration(narrow)) & ~isnan(ndepth(narrow)) );
          [r1,p1] = corr(nduration(narrow(bz)),ndepth(narrow(bz)),'type','Spearman');
          bz = find( ~isnan(nduration(broad)) & ~isnan(ndepth(broad)) );
          [r2,p2] = corr(nduration(broad(bz)),ndepth(broad(bz)),'type','Spearman');
          total = [narrow ; broad];
          tz = find( ~isnan(nduration(total)) & ~isnan(ndepth(total)) );
          [rt,pt] = corr(nduration(total(tz)),ndepth(total(tz)),'type','Spearman');
          title(sprintf('Nw:%4.2f(%5.3f);Bd:%4.2f(%5.3f)',r1,p1,r2,p2));
          %******
          disp('   ');
          disp(sprintf('%s : Depth correlations',topname));
          disp(sprintf('Nw: %6.4f(%8.6f); Bd:%6.4f(%8.6f); Tot:%6.4f(%8.6f)',r1,p1,r2,p2,rt,pt));
          %******
        else
          title(sprintf('Narrow(%d);Broad(%d)',length(narrow),length(broad)));
        end
        subplot('Position',[xx 0.1 0.20 0.20]);
        if ~correlated
           yx = hist(nduration,vx);
           bar(vx,yx,1,'Linestyle','none','FaceColor',[0,0,0],'FaceAlpha',0.3); hold on;  
        else
          nyx = hist(nduration(narrow),vx);
         % nyx = nyx / sum(nyx);
          byx = hist(nduration(broad),vx);
         % byx = byx / sum(byx);
          % plot(vx,nyx,'k.-','Color',ncolo,'Linewidth',2); hold on;
          %  plot(vx,byx,'k.-','Color',bcolo,'Linewidth',2);
          % bar(vx,nyx,1,'Linestyle','none','FaceColor',ncolo,'FaceAlpha',0.4); hold on;
            hold on;
            bar(vx,byx,1,'Linestyle','none','FaceColor',bcolo,'FaceAlpha',0.2); hold on;
          % plot([xmen,xmen],[V(3),V(4)],'k--');
        end
%         axis tight;
%         V = axis;
%         axis([xmin xmax V(3) (1.1*V(4))]);
%         V = axis;
        axis([-0.55 0.55 0 45]);
        V = [-0.55 0.55 0 45];
        xmin = -0.55;
        xmax = 0.55;
        if (correlated ~= 2) % && (correlated ~= 0)
            if (~correlated)
                xmen = 10.0;  % narrow vs broad thresh
            end
%            plot([xmen,xmen],[V(3),(1.1*V(4))],'k-');
             plot([0,0],[0 45],'k-');
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
        p = ranksum(nduration(narrow),nduration(broad));
        p1 = signrank(nduration(narrow));
        p2 = signrank(nduration(broad));
        pa = signrank(nduration);
        if (correlated)
            title(sprintf('N:%5.2f,B:%5.2f(p=%4.3f)A:%5.2f(p=%4.3f)',...
                   mod1,mod2,p,moda,pa));
            text(xmin+0.65*(xmax-xmin),0.90*V(4),sprintf('N:%5.2f(p=%4.3f)',mod1,p1),'Fontsize',8,'Color',ncolo);
            text(xmin+0.65*(xmax-xmin),0.80*V(4),sprintf('B:%5.2f(p=%4.3f)',mod2,p2),'Fontsize',8,'Color',bcolo);   
        
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
        xlabel(topname);
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
           %******************** 
        end
        %****************
        
return;



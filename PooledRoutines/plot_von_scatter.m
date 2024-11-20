function plot_von_scatter(AUC,LockInt,WithLog,WithStd,Range,color)

   %going to need three scatters, one for each parameter
   
    ACounts = zeros(1,4);
    cmap = [[0.4,0,0.4];[0.7,0.3,0.7];[0,0.4,0];[0.3,0.7,0.3];[0.99,0.73,0.01];[0.96,0.815,0.42]];
    for k = 1:size(AUC,1)
        colo = cmap(color,:);
        %****** color code the dots      % for now, plot as one color
        if (AUC(k,1) == 1)
            if (AUC(k,2) == 1)
               % colo = cmap(1,:);
                ACounts(1) = ACounts(1)+1;
            else
               % colo = cmap(2,:);
                ACounts(2) = ACounts(2)+1;
            end
        else
            if (AUC(k,1) == 3)
               if (AUC(k,2) == 1)
                %   colo = cmap(3,:);
                   ACounts(3) = ACounts(3)+1;
               else
                 %  colo = cmap(4,:);
                   ACounts(4) = ACounts(4)+1;
               end
            else
               if (AUC(k,2) == 1)
                  % colo = cmap(5,:);
                   ACounts(1) = ACounts(1)+1;
               else
                  % colo = cmap(6,:);
                   ACounts(2) = ACounts(2)+1;
               end
            end
        end
        %*******
        if (WithLog)
           h = loglog(AUC(k,5),AUC(k,3),'k.'); hold on;
           set(h,'Markersize',12);
           set(h,'Color',colo);
           if (WithStd)
             h2 = loglog([(AUC(k,5)-AUC(k,6)),(AUC(k,5)+AUC(k,6))],[AUC(k,3),AUC(k,3)],'k-');
             set(h2,'Color',colo);
             h3 = loglog([AUC(k,5),AUC(k,5)],[(AUC(k,3)-AUC(k,4)),(AUC(k,3)+AUC(k,4))],'k-');
             set(h3,'Color',colo);
           end
        else
           h = plot(AUC(k,5),AUC(k,3),'k.'); hold on;
           set(h,'Markersize',12);
           set(h,'Color',colo);
           if (WithStd)
             h2 = plot([(AUC(k,5)-AUC(k,6)),(AUC(k,5)+AUC(k,6))],[AUC(k,3),AUC(k,3)],'k-');
             set(h2,'Color',colo);
             h3 = plot([AUC(k,5),AUC(k,5)],[(AUC(k,3)-AUC(k,4)),(AUC(k,3)+AUC(k,4))],'k-');
             set(h3,'Color',colo);
           end
        end
        %*********
    end    
    axis tight;
    V = axis;
    if isempty(Range)
      maxa = max(V(2),V(4));
      maxa = maxa*1.25;
      mina = min(V(1),V(3))*0.5;
    else
       maxa = Range(2);
       mina = Range(1);
    end
    axis([mina maxa mina maxa]);
    if (WithLog)
       loglog([mina,maxa],[mina,maxa],'k-'); % unity line
    else
       plot([mina,maxa],[mina,maxa],'k-'); % unity line
    end
    %****** comp basic stats
    ai1 = (AUC(:,3) - AUC(:,5)) ./ (AUC(:,3)+AUC(:,5));
    [p1,h1] = signrank(ai1);
    medo = nanmedian(ai1);
    modo = (1+medo)/(1-medo);
    title(sprintf('Int(%4.2f,%4.2f) Mod:%5.3f (p=%6.4f)',...
                  LockInt(1),LockInt(2),modo,p1));
return;



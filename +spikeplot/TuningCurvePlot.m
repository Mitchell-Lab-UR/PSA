function TuningCurvePlot(OriTune,colo,symb,baroff,sidesmooth)
%*** function TuningCurvePlot(OriTune,colo,symb)
%*  OriTune - struct with mean response at each orientation
%*          - the orientation values
%*          - the fit curve to it, mu and sem
%   Colo   - desired color of the plot 
%   Symb   - symbol, usually '.' or 'o'
%   baroff - offset on x axis the error bars
%   sidesmooth - if 1, smooth adjacent, else raw
%
%*  Made this a function because later we can design it to center curve    

     %***********
     if isempty(OriTune)
         return;
     end
     if isempty(OriTune.Atune)
         return;
     end
     if isempty(OriTune.Atune.mu)
         return;
     end
     %***** recenter plots based on the peak ... write for later
     h = plot(OriTune.OriVals,OriTune.Atune.mu,'k-'); hold on;
     set(h,'Linewidth',2);
     set(h,'Color',colo);
     h = plot(OriTune.OriVals,(OriTune.Atune.mu + (2 * OriTune.Atune.sem)),'b-');
     fcolo = (0.75 * [1,1,1]) + (0.25 * colo);
     set(h,'Color',fcolo);
     h = plot(OriTune.OriVals,(OriTune.Atune.mu - (2 * OriTune.Atune.sem)),'b-');
     set(h,'Color',fcolo);
     %********
     if isfield(OriTune,'AAtune') && ~isempty(OriTune.AAtune)
         
       h = plot(OriTune.OriVals,OriTune.AAtune.mu,'k--'); hold on;
       set(h,'Linewidth',2);
       set(h,'Color',colo);
       h = plot(OriTune.OriVals,(OriTune.AAtune.mu + (2 * OriTune.AAtune.sem)),'b--');
       fcolo = (0.75 * [1,1,1]) + (0.25 * colo);
       set(h,'Color',fcolo);
       h = plot(OriTune.OriVals,(OriTune.AAtune.mu - (2 * OriTune.AAtune.sem)),'b--');
       set(h,'Color',fcolo);
     end
     
     %***** smooth raw values by next neighbor on circle ***
     zmu = OriTune.Mu;
     zsu = OriTune.Sem;
     if (sidesmooth == 1)
       LK = length(OriTune.OriVals);
       zzmu = zmu;
       zzsu = zsu;
       for k = 1:LK
         if (k==1)
            kset = [1,2,LK]; 
         else
             if (k == LK)
                kset = [(LK-1),LK,1];     
             else
                kset = [(k-1),k,(k+1)]; 
             end
         end
         zzmu(k) = mean(zmu(kset));
         zzsu(k) = mean(zsu(kset))/sqrt(3);
       end
       zmu = zzmu;
       zsu = zzsu;
     end
     %*********
     for k = 1:length(OriTune.OriVals)
       h = plot(OriTune.OriVals(k)+baroff,zmu(k),['k',symb]); hold on;
       set(h,'Markersize',10);
       set(h,'Color',colo*0.8);
       h = plot([OriTune.OriVals(k)+baroff,OriTune.OriVals(k)+baroff],...
                                   [zmu(k) + (2 * zsu(k)),...
                                    zmu(k) - (2 * zsu(k))],'k-'); hold on;
       set(h,'Color',colo*0.8);
     end
     %*****
     
  return;
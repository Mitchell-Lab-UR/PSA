%% Load Pooled Data
%Should live in the same path that you're running this script from 
%ie. C:\Users\amybu\Box\MTC\Results
     load('Pooled_VonMises_Results2_MTC_depth_Preload_191124.mat')
%%%

if (1)
    %zzt = find( (RATE_Anim > 0) & (R2List == 1));  % only good R2 fits
    Mcolo = [0,0.6,0.6];
    Ecolo = [0.3,0.7,0.3];
    h10 = figure(10); 
    set(h10,'Position',[100 100 1000 400]);
                  %subplot('Position',[0.1 0.15 0.35 0.70]); hold off;
    subplot (1,2,2)
    x = find (R2List == 1 & VONKAP{1}(:,1) == 1 &  VONKAP{1}(:,2) == 1); %Good Von Mises and MT and Sig.  
    x2 = find (R2List == 1 & VONKAP{1}(:,1) == 2 & VONKAP{1}(:,2) == 1); %Good Von Mises and MTC and Sig.  
    med = nanmedian(VONKAP{1}(x,5));
    med2 = nanmedian(VONKAP{1}(x2,5));

    %Plotting Histograms of Tuning Half Width
    histogram (VONKAP{1}(x,5), 14, 'Normalization', 'probability', 'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
    histogram (VONKAP{1}(x2,5), 14, 'Normalization', 'probability', 'FaceColor',[0.1 0.85 0.1], 'FaceAlpha', 0.5,  'EdgeColor', 'none'); hold on;
    xline (med,'Color',[0.8500 0.3250 0.0980] , 'LineWidth', 4)
    xline (med2, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 4)
    xlabel ('Tuning Half Width (deg)', 'FontSize', 18)
    ylabel ('Percentage of Units', 'FontSize', 18)
    title('Motion Tuning Curve Width','FontSize', 14)
     
     %Stats
     ai1 = VONKAP{1}(x,5)-90; %shift distribition around zero for wilcoxon signed rank test
     ai2 = VONKAP{1}(x2,5)-90;
     p1 = signrank(ai1);
     p2 = signrank(ai2);
     p = ranksum(ai1,ai2);
   % p = ranksum(VONKAP{1}(x,5), VONKAP{1}(x2,5));
     disp(sprintf('MTC:%5.2f(p=%5.3f);MT:%5.2f(p=%5.3f)',med,p1,med2,p2));
     disp(sprintf('Diff between areas p = %10.8f',p));
      
    if (0) %double check
     scatter (VONKAP{1}(x,5), VONKAP{1}(x,3) , '.')
    end

%%   DSI Histogram

    %Require Good Von Mises Fit
    s = find (R2List == 1 & DSI{1} > 0 & VONKAP{1}(:,1) == 1); %DSI greater than zero and MT
    s2 = find (R2List == 1 & DSI{1} > 0 & VONKAP{1}(:,2) == 1);
    
    %Just require DSI greater than zero (non-tuned units set at zero
    %s = find (DSI{1} > 0 & VONKAP{1}(:,1) == 1); %DSI greater than zero and MT
    %s2 = find (DSI{1} > 0 & VONKAP{1}(:,2) == 1);
    
    mes = nanmedian(DSI{1}(s));
    mes2 = nanmedian(DSI{1}(s2));
    
    %Plotting Histograms of Motion Tuning
    subplot (1,2,1)
    histogram (DSI{1}(s), 14, 'Normalization', 'probability', 'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
    histogram (DSI{1}(s2), 14, 'Normalization', 'probability', 'FaceColor',[0.1 0.85 0.1], 'FaceAlpha', 0.5,  'EdgeColor', 'none'); hold on;
    xline (mes,'Color',[0.8500 0.3250 0.0980] , 'LineWidth', 4)
    xline (mes2, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 4)
    xlabel ('DSI', 'FontSize', 18)
    ylabel ('Percentage of Units', 'FontSize', 18)
    title('Motion Tuning Strength','FontSize', 14)
    
    %Stats
    d1 = DSI{1}(s); %shift distribition around zero for wilcoxon signed rank test
    d2 = DSI{1}(s2);
    p1 = signrank(d1);
    p2 = signrank(d2);
    p = ranksum(d1,d2);
    disp(sprintf('MTC:%5.2f(p=%5.3f);MT:%5.2f(p=%5.3f)',mes,p1,mes2,p2));
    disp(sprintf('Diff between areas p = %10.8f',p));
%%%

end
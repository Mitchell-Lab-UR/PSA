function [vec] = Compute_fstat_fano_rate_matched(xx1,yy1,xx2,yy2,CWin)

    fano = NaN;
    
    %***** compute the f-stat, between var by within var
    mu1 = [];
    va1 = [];
    xvals = unique(xx1);
    for ix = 1:length(xvals)
        zz = find( xx1 == xvals(ix));
        if (length(zz) > 1)
          mu1 = [mu ; mean(yy1(zz))];  
          va1 = [va ; var(yy1(zz))];   
        end
    end
    mu1 = mu1 * CWin;   % back to counts
    va1 = va1 * (CWin^2); % back to counts
    %******* then do the second condition
    mu2 = [];
    va2 = [];
    xvals = unique(xx2);
    for ix = 1:length(xvals)
        zz = find( xx2 == xvals(ix));
        if (length(zz) > 1)
          mu2 = [mu2 ; mean(yy2(zz))];  
          va2 = [va2 ; var(yy2(zz))];   
        end
    end
    mu2 = mu2 * CWin;   % back to counts
    va2 = va2 * (CWin^2); % back to counts
    %******** now, try to find a matching distribution in mu vs mu2
    if (1)
        figure(10);
        plot(mu1,va1,'ro'); hold on;
        plot(mu2,va2,'bo');
        xlabel('mean');
        ylabel('variance');
        input('check'); 
    end
    zz = 1:length(mu1);
    %**********
    fano1 = mean( (va1(zz) ./ mu1(zz)) );  % more sensitive measure
    fano2 = mean( (va2(zz) ./ mu2(zz)) );  % more sensitive measure
     
    vec = [fano1,fano2];
    
return;
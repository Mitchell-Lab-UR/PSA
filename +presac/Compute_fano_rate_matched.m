function [vec] = Compute_fstat_fano_rate_matched(xx1,yy1,xx2,yy2,CWin)

    fano = NaN;
    
    %***** compute the f-stat, between var by within var
    mu1 = [];
    va1 = [];
    xvals = unique(xx1);
    for ix = 1:length(xvals)
        zz = find( xx1 == xvals(ix));
        if (length(zz) > 1)
          mu1 = [mu1 ; mean(yy1(zz))];  
          va1 = [va1 ; var(yy1(zz))];   
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
    %***** matching process here
    %** go through attended list, find closest matching point and
    %** accept if rate is withing 5 percent ... continue from there
    zz1 = [];
    zz2 = [];
    zgone = 1:length(mu2);
    for zi = 1:length(mu1)
        amu = mu1(zi);
        dist = abs( mu2(zgone) - amu);
        z = find( dist == min(dist));
        zj = zgone(z(1));
        bmu = mu2(zj);
        if ((amu+bmu) == 0)
            continue;
        else
            pmod = (amu-bmu)/(amu+bmu);
            if (abs(pmod) < 0.03)  % within 5 percent mod either up or down
                zz1 = [zz1 zi];
                zz2 = [zz2 zj];
                zg = find( zgone == zj);
                zgone = zgone([1:(zg-1),(zg+1):end]);  % omit from further tests
            end
        end
    end
    if (length(zz1) < 4)  % require enough points were matched, otherwise forget it
        fano1 = NaN;
        fano2 = NaN;
        fanoS1 = NaN;
        fanoS2 = NaN;
    else
        fano1 = mean( (va1(zz1) ./ mu1(zz1)) );  % more sensitive measure
        fano2 = mean( (va2(zz2) ./ mu2(zz2)) );  % more sensitive measure
        fanoS1 = regress(va1(zz1),mu1(zz1));
        fanoS2 = regress(va1(zz1),mu1(zz1));
    end
    %******** now, try to find a matching distribution in mu vs mu2
    if (0)
        figure(10); hold off;
        plot(mu1,va1,'ro'); hold on;
        plot(mu1(zz1),va1(zz1),'r.');
        plot(mu2,va2,'bo');
        plot(mu2(zz2),va2(zz2),'b.');
        xlabel('mean');
        ylabel('variance');
        [mean(mu1),mean(mu2)]
        [mean(mu1(zz1)),mean(mu2(zz2))]
        input('check'); 
    end
    %**********
     
    vec = [fano1,fano2,fanoS1,fanoS2];
    
return;
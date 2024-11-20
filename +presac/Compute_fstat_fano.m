function [vec] = Compute_fstat_fano(xx,yy,CWin)

    finfo = NaN;
    fano = NaN;
    
    %***** compute the f-stat, between var by within var
    mu = [];
    va = [];
    xvals = unique(xx);
    for ix = 1:length(xvals)
        zz = find( xx == xvals(ix));
        if (length(zz) > 1)
          mu = [mu ; mean(yy(zz))];  
          va = [va ; var(yy(zz))];   
        end
    end
    finfo = sqrt(var(mu)/mean(va));   % between sigma divided by within 
                                      % reduces to d-prime for two values + Gauss
    %***** back to counts before doing Fano Factor
    mu = mu * CWin;   % back to counts
    va = va * (CWin^2); % back to counts
    %*********** fit slope line for mean vs var ... fano-factor
    fano = regress(va,mu);
    zz = find( mu > 0);
    fano2 = mean( (va(zz) ./ mu(zz)) );  % more sensitive measure
    
    if (0) % check plot for sanity sake
      y = fano * mu;
      figure(100);
      plot(mu,va,'ko'); hold on;
      plot(mu,y,'b-');
      axis square;
      [finfo,fano]
      input('check');
    end
    
    vec = [finfo,fano,fano2];
    
return;
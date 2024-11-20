function [AttWid,IgnWid] = compute_half_width(mu1,mu2)
                
        ka = mu1(3);  % kappa parameter
        if (ka > 0)
           AttWid = acos(log(.5 + .5*exp(2*ka))./ka -1) * (180/pi);
        else
           AttWid = 180 - (acos(log(.5 + .5*exp(2*abs(ka)))./abs(ka) -1) * (180/pi));    
        end
        kb = mu2(3);  % kappa parameter
        if (kb > 0)
           IgnWid = acos(log(.5 + .5*exp(2*kb))./kb -1) * (180/pi);
        else
           IgnWid = 180 - (acos(log(.5 + .5*exp(2*abs(kb)))./abs(kb) -1) * (180/pi));    
        end
        
return;
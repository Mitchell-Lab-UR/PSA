function sy = circsmooth(y,np)
   sy = y;
   np = floor(np);  % must be integer
   if (np == 0)
       return;
   else
     for ii = 1:length(y)
        ia = ii-np;
        ib = ii+np;
        sumo = [];
        for jj = ia:ib
            it = jj;
            if (it < 1)
                it = it + length(y);
            end
            if (it > length(y))
                it = it - length(y);
            end
            sumo = [sumo y(it)];
        end
        sy(ii) = nanmean(sumo);
     end
   end
return;

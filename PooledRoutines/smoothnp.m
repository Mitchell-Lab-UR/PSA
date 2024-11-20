function sy = smoothnp(y,np)
   sy = y;
   np = floor(np);  % must be integer
   if (np == 0)
       return;
   else
     for ii = 1:length(y)
         ia = max(1,(ii-np));
         ib = min(length(y),(ii+np));
         sy(ii) = nanmean(y(ia:ib));
     end
   end
return;

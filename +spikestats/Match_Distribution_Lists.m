function BigList = Match_Distribution_Lists(AttList,IgnList,RTList,BinSize)
%  function BigList = Match_Distribution_Lists(AttList,IgnList,RTList,BinSize)
%   inputs:
%      AttList - list of trials for attention condition
%      IgnList - list of trials for unattended condition
%      RTList - list of RTs per trial
%      BinSize - size of bins for matching (if RTs, 0.01 secs for example)
%
%   outputs:   Cell(1,2) - list of trials matched in parameter

   attrts = RTList(AttList);
   ignrts = RTList(IgnList);
   mino = min(min(attrts),min(ignrts));
   cmino = floor((mino/BinSize)) * BinSize;
   maxo = max(max(attrts),max(ignrts));
   cmaxo = ceil((maxo/BinSize)) * BinSize;
   vx = cmino:BinSize:cmaxo;   
   hx = (cmino+(BinSize/2)):BinSize:(cmaxo-(BinSize/2));
   %***** make rate matched lists at random draw
   amat = [];
   imat = [];
   for k = 1:(length(vx)-1)
       aa = vx(k);
       bb = vx(k+1);
       za = find( (attrts >= aa) & (attrts < bb));
       zi = find( (ignrts >= aa) & (ignrts < bb));
       if (length(za) > length(zi))
           rr = randperm(length(za));
           za = za(rr(1:length(zi)));
       else
           if (length(za) < length(zi))
               rr = randperm(length(zi));
               zi = zi(rr(1:length(za)));
           end
       end
       if ~isempty(za)
           amat = [amat ; AttList(za)];
       end
       if ~isempty(zi)
           imat = [imat ; IgnList(zi)];
       end
   end
   BigList = cell(1,2);
   BigList{1} = amat;
   BigList{2} = imat;
   
   %******** then plot the new distributions
   if (0)
      atty = hist(attrts,hx);
      itty = hist(ignrts,hx);
      %****** recompute dists
      attrts2 = RTList(amat);
      ignrts2 = RTList(imat);
      atty2 = hist(attrts2,hx);
      itty2 = hist(ignrts2,hx);
      %******* 
      figure(1);
      subplot(2,1,1); 
      plot(hx,atty,'r.-'); hold on;
      plot(hx,itty,'bo:');
      xlabel('rt (secs)');
      ylabel('unmatched counts');
      subplot(2,1,2); 
      plot(hx,atty2,'r.-'); hold on;
      plot(hx,itty2,'bo:');
      xlabel('rt (secs)');
      ylabel('matched counts');
      input('stop');
   end
   
end



function pairedcorr = compute_noise_corr(CorrOri,CorrSpk,CorrSig,CorrTyp)
   pairedcorr = [];
   M = size(CorrOri{1,1},1);  % number of units in cluster
   %***** first, along each row remove the mean spike count per motion
   %***** motion direction to eliminate any signal in the correlations
   for zk = 1:2
       ovals = unique(CorrOri{1,zk});
       for k = 1:M
         for oo = 1:length(ovals)
             zz = find( CorrOri{1,zk}(k,:) == ovals(oo));
             if (length(zz) > 1)
               uu = nanmean( CorrSpk{1,zk}(k,zz) );
               su = nanstd( CorrSpk{1,zk}(k,zz) );
               if (su > 0)
                 CorrSpk{1,zk}(k,zz) = (CorrSpk{1,zk}(k,zz) - uu) /su;  % z score
               else
                 CorrSpk{1,zk}(k,zz) = 0;  % if su == 0, then all values = 0   
               end
             else
                 CorrSpk{1,zk}(k,zz) = NaN;  % throw out if not defined      
             end
         end
       end
   end
   %********** then 
   for k = 1:M
       for j = (k+1):M
            ncorr = [];
            nrange = [];
            for zk = 1:2  % 1 is towards and 2 is away
                zz = find( ~isnan(CorrSpk{1,zk}(k,:)) & ~isnan(CorrSpk{1,zk}(j,:)) );
                if (length(zz) > 2)
                  [r,p,rl,ru] = corrcoef(CorrSpk{1,zk}(k,zz)',CorrSpk{1,zk}(j,zz)');
                  ncorr = [ncorr r(1,2)];
                  nrange = [nrange (ru(1,2)-rl(1,2))];
                  % ncorr = [ncorr corr(CorrSpk{1,zk}(k,zz)',CorrSpk{1,zk}(j,zz)')]; %,'type','Spearman')];
                else
                  ncorr = [];
                  break;
                end
            end
            nsig = corr(CorrSig{1,1}(k,:)',CorrSig{1,1}(j,:)'); %,'type','Spearman');
            %*****
            if ~isempty(ncorr)
                if (ncorr(1) == 1) && (ncorr(2) == 1) % impossible, means a unit was double counted?
                   disp(sprintf('check here, there is a problem in unit double count'));
                   % input('check');
                else
                   pairedcorr = [pairedcorr ; [ncorr nsig nrange CorrTyp{1,1}(k) CorrTyp{1,1}(j)]];
                end
            end
       end
   end
   %************
return;

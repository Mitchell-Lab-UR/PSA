function h = PlotTickRaster(Rast,TrialSort,Voffset)
%*** function PlotTickRaster(Rast,TrialSort)
%**
%** inputs:  Rast:  Nx2 where column 1 is times and column 2 is trials
%**          TrialSort:  if included, a list of variables per trial over
%**                 which to sort the vertical scale of rasters
%**          Voffset:  draw rasters at some vertical offset
%** outputs:   h - vectors of tick marks (change color, or size)
%**
%**   uses trick from Jack to plot tick marks better
%**   plots in the current subplot vertical lines as spike ticks
%**   returns the handle to all the tick marks (if change color etc..)

%****** if too few spikes, don't even bother
if isempty(Rast) || (size(Rast,1) < 10)
    h = plot([1,1],[1,2],'k.-');
    return;
end

%** use a sort list to reorder trials vertically
if ~isempty(TrialSort)
  trialids = 1:length(TrialSort);
  tuningids = TrialSort;
  if (size(tuningids,2) < size(tuningids,1) )
      tuningids = tuningids';
  end
  sorto = [trialids' tuningids'];
  ysort = sortrows(sorto,2);
  %***************
  for k = 1:length(Rast)
      tr = Rast(k,2);  % trial tag
      z = find( ysort(:,1) == tr);
      Rast(k,2) = z(1);
  end
  %**************
end

%*********
ix = Rast(:,1)';
iy = Rast(:,2)';

n = numel(ix);
height = 1;

x = [ix; ix; nan(1,n)];
y = [iy; iy+height; nan(1, n)];
    
h = plot(x(:), (y(:)+Voffset), 'k-');



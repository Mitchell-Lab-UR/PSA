function plot_attention_traces(TT,UU,SU,LockWin,LockInt,Event)  
%** function plot_attention_traces(TT,UU,SU,LockWin,Event)  
%**   
%**    routine to plot PSTH, or other time traces (Att vs Ignored)
%**    assumes two conditions, 1 is att, 2 is ignored for TT, UU, SU as
%**             cell structs 

 for zk = 1:2    %could do 1:3 if want to see less sampled ignored loc
     if (zk == 1)
         colo = [1,0,0]; % att in red
     else
         colo = [0,0,1]; % ign1 in blue
     end
     tt = TT{zk};
     uu = UU{zk};
     su = SU{zk};
     %********* now plot results
     h2 = plot(tt,uu,'k-');  hold on;
     set(h2,'Color',colo*0.8);
     set(h2,'Linewidth',2);
     tta = [tt fliplr(tt)];
     yya = [(uu+(2*su)) fliplr((uu-(2*su)))];
     fill(tta,yya,colo,'FaceAlpha',0.3,'Linestyle','none');
   end
   axis tight;
   V = axis;
   VM = V(4);
   axis([LockWin(1) LockWin(2) 0 VM]);
   % plot([0,0],[0,VM],'k-');
   if ~isempty(LockInt) % use gray fill to mark time window
      aa = [LockInt(1),LockInt(1),LockInt(2)-0.001,LockInt(2)-0.001];
      bb = [0,VM,VM,0];
      fill(aa,bb,[0.5,0.5,0.5],'FaceAlpha',0.3,'Linestyle','none');
   end
   if isempty(Event)
       xlabel('Stimulus Rank Order','Fontsize',15);
   else
     if (Event == 1)
       TiName = 'Stimulus Onset';
     end
     if (Event == 2)
       TiName = 'Saccade Onset';
     end
     if (Event == 3)
       TiName = 'Saccade Offset';
     end
      xlabel(['Time from ',TiName,'(s)'], 'FontSize', 15);
   end
   set(gca,'FontSize',15);
 return;
   
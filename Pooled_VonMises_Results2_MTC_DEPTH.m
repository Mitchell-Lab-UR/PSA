function Pooled_VonMises_Results2()
%*** grabs results from Info directory, goes unit by unit and
%*** applies selection criteria, creates pooled results


%***
if (0) % Jude's Laptop
   StorePath = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Info_Sprout';  % MTC Laminar Sprout
   StorePath2 = 'C:\Users\jmitchell\Box Sync\FovealTrans\MTC\Info_MTC';  % MTC Laminar Milo
else
    StorePath = 'C:\Users\amybu\Box\MTC\Info_MTC';
   %StorePath = 'C:\Users\amybu\Box\MTC\Info_MTC';  % MTC Laminar
   StorePath2 = 'C:\Users\amybu\Box\MTC\Info_MT';  % MT Laminar
    
end

if (1)%*** Load Depth information from Gabe's CSDs


if (1)
        
   % *****************MT ONLY
         if(1) %MT
    
          LaminarNames1 = {'Milo_301121','Milo_021221','Milo_161221','Milo_201221',...
                          'Milo_040122','Milo_060122','Milo_180122','Milo_280122',...
                          'Milo_010222','Milo_080222','Milo_250222','Milo_010322',...
                          'Milo_180322','Milo_050422',...  %add Milo MT
                          'Sprout_080922','Sprout_140922','Sprout_141222',...
                          'Sprout_161222','Sprout_281222','Sprout_040123',...
                          'Sprout_250123', 'Sprout_270123',...
                          'Sprout_010223','Sprout_030223','Sprout_100223',... %Sprout MT
                          };


          LaminarNames1 = string(LaminarNames1)';
          Laminar1Depths_Shank1 = [ 525	945 NaN NaN NaN NaN NaN NaN NaN 550 700 630 NaN 400 ...
                                     400 420 570 560 665 665 700 700 700 580 525]; % (MT SPROUT) This is the bottom of the input layer 

          Laminar1Depths_Shank2 = [525	875 NaN NaN NaN NaN NaN NaN NaN 620 700 630 NaN 400 ...
                                    420	420 570	560	665	735 735	735	700	580	630]; %MT SPROUT
          InputWidth = 300; 
         else
            % 
            % %**************Foveal MT and MTC MIX
            %   LaminarNames = {'Milo_190221','Milo_190421','Milo_230421','Milo_260421',...
            %                   'Milo_300421','Milo_030521','Milo_120521','Milo_280521',...
            %                   'Milo_030621','Milo_210721','Milo_220721','Milo_120821',...
            %                   'Ellie_090120','Ellie_271219','Ellie_180120','Ellie_190120',...
            %                   'Ellie_200120','Ellie_210120', 'Ellie_230120', 'Ellie_270120'...
            %                   'Ellie_250220', 'Ellie_270220'};
            %   LaminarNames = string(LaminarNames)';
            %   LaminarDepths_Shank1 = [600 805 800 775 NaN 520 665 980 770 735 650 735 665 630 725 770 770 805 800 780 280 NaN 800]; % This is the bottom of the input layer 
            %   LaminarDepths_Shank2 = [600 735 805 725 NaN 595 700 945 700 735 650 700 665 525 725 665 665 820 875 875 385 NaN 800];
            %   InputWidth = 300; 
            %   %******* 
         end
  
    %else 
        
        if (1) %MTC
            LaminarNames2 = {'Milo_300921', 'Milo_051021', 'Milo_071021','Milo_121021',...
                       'Milo_141021', 'Milo_211021', 'Milo_091121', 'Milo_111121',...
                       'Milo_141221', 'Milo_201221', 'Milo_180222', ...
                       'Milo_220222', 'Milo_220322', 'Milo_290322', ...
                       'Milo_010422', 'Milo_080422', 'Milo_150422','Milo_190422',... %Milo MTC
                      'Sprout_300922','Sprout_051022','Sprout_121022',...
                      'Sprout_141022','Sprout_181022','Sprout_201022','Sprout_251022',...
                      'Sprout_031122','Sprout_071222', 'Sprout_091222',...
                      'Sprout_130123','Sprout_180123','Sprout_010323',... %Sprout MTC
                      };

          LaminarNames2 = string(LaminarNames2)';
          Laminar2Depths_Shank1 = [ 595	665	525	630	700	NaN	630 NaN 945 NaN 175 525 280	665	280 210 NaN 245 ...
                                   490	700	385	630	595	525	630	560	490	455 525	630 NaN]; % (MTC SPROUT) This is the bottom of the input layer 

          Laminar2Depths_Shank2 = [595	700	525	595	665	NaN	630 NaN 945 NaN 290 490 280	665	300 210 NaN 350 ...
                                   560	700	420	630	630	560	665	595	525	500 595	665 NaN]; % (MTC SPROUT)
          InputWidth = 300; 
          
        else
              % %*********************SHANNA OLD MT AND MTC MIX
              % %**************SHANNA OLD MT AND MTC MIX
              % LaminarNames = {'Milo_150221','Milo_220221','Milo_170321','Milo_180321',...
              %                 'Milo_070421','Milo_120421','Milo_160421','Milo_210421',...
              %                 'Milo_100521','Milo_140521','Milo_080721','Milo_130721',...
              %                 'Milo_150721','Milo_280721','Milo_300721','Milo_170821',...
              %                 'Milo_240821','Milo_260821', 'Milo_160921'};
              % LaminarNames = string(LaminarNames)';
              % LaminarDepths_Shank1 = [525 735 840 875 665 840 630 735 805 665 770 595 875 630 980 805 735 525 770]; % This is the bottom of the input layer 
              % LaminarDepths_Shank2 = [NaN NaN 770 840 700 875 595 700 805 630 805 595 840 875 945 805 700 560 770];
              % InputWidth = 300; 
              % %******* 
        end
    end

LaminarNames = vertcat(LaminarNames1,LaminarNames2);
LaminarDepths_Shank1 = horzcat(Laminar1Depths_Shank1, Laminar2Depths_Shank1 );
LaminarDepths_Shank2 = horzcat(Laminar1Depths_Shank2, Laminar2Depths_Shank2 );

end
%***********************************
%%
if (0)  %%Load Depth information from CSDs (Copied from SacMod Script)
    
    
    if ~MTC
        
    
    
  LaminarNames = {'Milo_301121','Milo_021221','Milo_161221','Milo_201221',...
                  'Milo_040122','Milo_060122','Milo_180122','Milo_280122',...
                  'Milo_010222','Milo_080222','Milo_250222','Milo_010322',...
                  'Milo_180322','Milo_050422',...  %add Milo MT
                  'Sprout_080922','Sprout_140922','Sprout_141222',...
                  'Sprout_161222','Sprout_281222','Sprout_040123',...
                  'Sprout_250123', 'Sprout_270123',...
                  'Sprout_010223','Sprout_030223','Sprout_100223',... %Sprout MT
                  };

%                   %Milo MTC
%                    'Milo_300921', 'Milo_051021', 'Milo_071021','Milo_121021',...
%                    'Milo_141021', 'Milo_211021', 'Milo_091121', 'Milo_111121',...
%                    'Milo_141221', 'Milo_201221', 'Milo_180222', ...
%                    'Milo_220222', 'Milo_220322', 'Milo_290322', ...
%                    'Milo_010422', 'Milo_080422', 'Milo_150422','Milo_190422'...
                  
%                   %Sprout MTC
%                   'Sprout_300922','Sprout_051022','Sprout_121022',...
%                   'Sprout_141022','Sprout_181022','Sprout_201022','Sprout_251022',...
%                   'Sprout_031122','Sprout_071222', 'Sprout_091222',...
%                   'Sprout_130123','Sprout_180123','Sprout_010323'
                  
           
              
  LaminarNames = string(LaminarNames)';
  LaminarDepths_Shank1 = [ 525	945 NaN NaN NaN NaN NaN NaN NaN 550 700 630 NaN 400 ...
                             400 420 570 560 665 665 700 700 700 580 525]; % (MT SPROUT) This is the bottom of the input layer 
                         
  LaminarDepths_Shank2 = [525	875 NaN NaN NaN NaN NaN NaN NaN 620 700 630 NaN 400 ...
                            420	420 570	560	665	735 735	735	700	580	630]; %MT SPROUT
  InputWidth = 300; 
  
    else 
        LaminarNames = {'Milo_300921', 'Milo_051021', 'Milo_071021','Milo_121021',...
                   'Milo_141021', 'Milo_211021', 'Milo_091121', 'Milo_111121',...
                   'Milo_141221', 'Milo_201221', 'Milo_180222', ...
                   'Milo_220222', 'Milo_220322', 'Milo_290322', ...
                   'Milo_010422', 'Milo_080422', 'Milo_150422','Milo_190422',... %Milo MTC
                  'Sprout_300922','Sprout_051022','Sprout_121022',...
                  'Sprout_141022','Sprout_181022','Sprout_201022','Sprout_251022',...
                  'Sprout_031122','Sprout_071222', 'Sprout_091222',...
                  'Sprout_130123','Sprout_180123','Sprout_010323',... %Sprout MTC
                  };
              
  LaminarNames = string(LaminarNames)';
  LaminarDepths_Shank1 = [ 595	665	525	630	700	NaN	630 NaN 945 NaN 175 525 280	665	280 210 NaN 245 ...
                           490	700	385	630	595	525	630	560	490	455 525	630 NaN]; % (MTC SPROUT) This is the bottom of the input layer 
                       
  LaminarDepths_Shank2 = [595	700	525	595	665	NaN	630 NaN 945 NaN 290 490 280	665	300 210 NaN 350 ...
                           560	700	420	630	630	560	665	595	525	500 595	665 NaN]; % (MTC SPROUT)
  InputWidth = 300; 
              
                   
                   
        
    end
  %******* 
  
end
%%
%*******
InfoNames = dir(StorePath);
InfoNames2 = dir(StorePath2);
N1 = size(InfoNames,1)-2;
N2 = size(InfoNames2,1)-2;
%N2=0;
%***
FntSize = 12;
%***** AUC columns: (animal) (iso) (att auc) (att sem) (ign att) (ign sem)
R2List = [];  % units with good enough VonMises to use in curve analyses
zR2List = [];  % for feature gain tests
AUC_tt = [];
AUC_att = [];
AUC_ign = [];
ATTAUC = [];
IGNAUC = [];
RATE_tt = cell(1,2);
RATE_att = cell(1,2);
RATE_ign = cell(1,2);
RATE_Anim = [];
VONBASE = cell(1,2);
VONGAIN = cell(1,2);
VONKAP = cell(1,2);
VON_Anim = [];
ZVONBASE = cell(1,2);
ZVONGAIN = cell(1,2);
ZVONKAP = cell(1,2);
%*****
VONFIT = cell(1,2);
VONFITSEM = cell(1,2);
ZVONFIT = cell(1,2);
ZVONFITSEM = cell(1,2);
DSI = cell (1,2);
%******
NDepth = [];
Layer = [];
CSD = [];
Duration = [];
%*******
FR = cell(1,2);
AUC = cell(1,2);
MI = cell(1,2);
MIC = cell(1,2);
MIC2 = cell(1,2);
FF = cell(1,2);
FS = cell(1,2);
%*** non-parametric tuning curves
NPSMOOTH = 1;  % adjacent bin smoothing
NPFIT = cell(1,2);
NPFEAT = cell(1,2);
NPAMP = [];
NPBASE = [];
NPWID = [];
NPMIN = [];  % find minimum mean rate per condition (criterion for incl?)
RNPAMP = [];
RNPBASE = [];
RNPWID = [];
LOCBIAS = [];

%***** going through files, look for pairs and accumulate counts
OldFile = [];
CorrSpk = cell(1,2);
CorrOri = cell(1,2);
CorrSig = cell(1,1);
CorrTyp = cell(1,1);
NCORR = [];  % will list Animal, Towards NCorr, Away NCorr (per pair)

%***** store all processed lists and info in environment so
%***** you don't need to rerun it each time
%Environment_Preload = 'Pooled_VonMises_Results2_MTC_ALL_Preload.mat';
Environment_Preload = 'Pooled_VonMises_Results2_MTC_depth_Preload_191124.mat';
showall = 1;
MINRATE = 1;  % if mean rate in unattended or attended sac onset is
              % less than this exclude (impossible to fit it well)
Rthresh = 0.5;   % Could we do better with log-likelihood ratio?
Wthreshmin = 9.5;
Wthreshmax = 11.5;
%**** cell counts
Miso = 0; Mmu = 0; Mtot = 0; Mtuned = 0;  Mtungsten = 0;  Mincrate = 0;
Eiso = 0; Emu = 0; Etot = 0; Etuned = 0;  Etungsten = 0;  Eincrate = 0;
Msingle = 0; Mmulti = 0;
Esingle = 0; Emulti = 0;
ExampleInfo = [];
ExampleName = 'MT_260821_U17.mat';  % low rate, best ever

 

%*****
centerwid = 0;
if (0)
%if exist(Environment_Preload)
    disp(sprintf('Loading precomputed lists from %s ...',Environment_Preload));
    load(Environment_Preload)
    disp('Loading completed');
    StorePath = 'C:\Users\amybu\Box\MTC\Info_MTC';  % MTC Laminar
    StorePath2 = 'C:\Users\amybu\Box\MTC\Info_MT';  % MT Laminar
    InfoNames = dir(StorePath);
    InfoNames2 = dir(StorePath2);
    N1 = size(InfoNames,1)-2;
    N2 = size(InfoNames2,1)-2;
    input 'check'

else

 
  %*********
  for k = 1:(N1+N2)  % grab files from two directories
      tungsten = 0; % all recordings are laminar here
      if (k <= N1)
          animal = 1;  % here tungsten can mean MTC folder
          fname = InfoNames(k+2).name;
          % continue;
          disp(sprintf('Processing tungsten unit %s',fname));
          load([StorePath,filesep,fname]);   % load Info struct
      else
            animal = 2;
            % continue;  % skip MT laminars for now
            fname = InfoNames2(k+2-N1).name;
            disp(sprintf('Processing laminar unit %s',fname));
            load([StorePath2,filesep,fname]);   % load Info struct
      end
      if (1)
          %*****
          if (Info.isolation(1) >= 1)
              iso = 1;  % single unit and cluster MUA
          else
              iso = 0;  % hash/noise
          end
          includener = 0;
          if isfield(Info,'Included')
              includener = 1;
          end
          
          %********** inc all for now
          % if (tungsten == 0)  % don't include hash from Kilosort of MT
          %   if (iso == 0)   % any from MTC OK, including MU
          %     includener = 0;
          %   end
          % end
          
          incrate = 0;
          if isfield(Info,'Included') && (Info.Included > 0)
             ratefit1 = Info.SacOn.OriTune{1}.Mu;  % orientation tuning
             ratefit2 = Info.SacOn.OriTune{2}.Mu;
             minrate = min([nanmean(ratefit1),nanmean(ratefit2)]);  % min rate in condition
             if (minrate < MINRATE)
                incrate = 0;
             else
                incrate = 1;
             end
          end
          
          %********* tally up stats ********
          if (animal == 1)
             Etot = Etot + 1;
             if isfield(Info,'Included') && (Info.Included > 0) && incrate
                 if (iso > 0)
                     Eiso = Eiso + 1;
                 else
                     Emu = Emu + 1;
                 end
                 if (Info.TunedUnit)
                     Etuned = Etuned + 1;
                     if (incrate)
                        Eincrate = Eincrate + 1;
                        if (iso > 0)
                            Esingle = Esingle + 1;
                        else
                            Emulti = Emulti + 1;
                        end 
                     end
                 end
             end
          else
             if (animal == 2)  
                Mtot = Mtot + 1;
                if isfield(Info,'Included') && (Info.Included > 0) && incrate
                   if (iso > 0)
                     Miso = Miso + 1;
                   else
                     Mmu = Mmu + 1;
                     Mtungsten = Mtungsten + tungsten;
                   end
                   if (Info.TunedUnit)
                     Mtuned = Mtuned + 1;
                     if (incrate)
                        Mincrate = Mincrate + 1;
                        if (iso > 0)
                            Msingle = Msingle + 1;
                        else
                            Mmulti = Mmulti + 1;
                        end 
                     end
                   end
                end 
             end
          end
      end
      %***********
      if ~isfield(Info,'Included')
          continue;
      end
      if (Info.Included > 0)
          %*****
          % if (Info.isolation(1) > 2)
          if (Info.isolation(1) >= 1)
              iso = 1;  % single unit and cluster MUA
          else
              iso = 0;  % hash/noise
          end
          %*******
          
          %**** minimum firing rate during stimulation conditions
          ratefit1 = Info.SacOn.OriTune{1}.Mu;  % orientation tuning
          ratefit2 = Info.SacOn.OriTune{2}.Mu;
          minrate = min([nanmean(ratefit1),nanmean(ratefit2)]);  % min rate in condition
          if (minrate < MINRATE)
              includener = 0;
          end

          %****** require significant DSI tuning *********
          %if ~Info.TunedUnit
          %    includener = 0;
          %end
          
          %**********
          if (~includener)
                continue;
          end

          %************** add laminar depth information here 
          if (fname(1) == 'M')
             x = find( (fname(1:11) == LaminarNames) == 1);   
          else
             x = find( (fname(1:13) == LaminarNames) == 1);
          end  
          if ~isempty(x) % include depth information
                  if (Info.shank == 1)
                      LaminarDepths = LaminarDepths_Shank1;
                  else
                      LaminarDepths = LaminarDepths_Shank2;
                  end
                  Info.CSD = LaminarDepths(x);
                  depth = Info.CSD-Info.depth;
                  if (depth < 0)
                      Info.layer = 3;  % deep
                  else
                      if (depth > InputWidth)
                          Info.layer = 1; % superficial
                      else
                          Info.layer = 2; % input layer
                      end
                  end
                  %*****
          else
             %disp(sprintf('Cant find layer for %s',fname));
             Info.layer = NaN;
             Info.CSD = NaN;
             Info.depth = NaN;
          end
          %**************
          
          if (animal == 2) %******* COMPUTE PAIRWISE NOISE CORRELATIONS  
            %**** if neuron included, then check for noise correlation
            tagname = Info.tagname;
            z = find(tagname == 'U');
            tagname = tagname(1:(z(1)-2));  % name before _U$
            %****** check if matches old file
            if ~isempty(OldFile) && ~strcmp(tagname,OldFile)  % new cluster set
                N = size(CorrOri{1,1},1);  % how many rows?
                if (N > 1)   % compute ncorrs and add to NCORR
                   pairedcorr = compute_noise_corr(CorrOri,CorrSpk,CorrSig,CorrTyp);
                   M = size(pairedcorr,1);
                   for am = 1:M
                       NCORR = [NCORR ; [animal pairedcorr(am,1:7)]];
                   end
                end
                %****** reset things
                CorrOri = cell(1,2);
                CorrSpk = cell(1,2);
                CorrSig = cell(1,1); % use tuning curves (signal corr)
                CorrTyp = cell(1,1);
                %**********
            end
            OldFile = tagname;
            %***** regardless add new spike counts to lists
            for ak = 1:2   %note, if not matched in counts will be an error
              % CorrOri{1,ak} = [CorrOri{1,ak} ; Info.AUC_SacOn.SpkOri{ak,1}'];
              % CorrSpk{1,ak} = [CorrSpk{1,ak} ; Info.AUC_SacOn.SpkCnt{ak,1}'];
              CorrOri{1,ak} = [CorrOri{1,ak} ; Info.AUC_SacOff.SpkOri{ak,1}'];
              CorrSpk{1,ak} = [CorrSpk{1,ak} ; Info.AUC_SacOff.SpkCnt{ak,1}'];
            end
            CorrSig{1,1} = [CorrSig{1,1} ; Info.OriTune.Mu];  % motion tuning
            %******** set category, layer x duration = 3 x 2 = 6
            nb = NaN;
            if (Info.duration < Wthreshmin)
                nb = 1;
            end
            if (Info.duration > Wthreshmin)
                nb = 2;
            end
            CorrTyp{1,1} = [CorrTyp{1,1} ; (((Info.layer-1)*2) + nb)];
            %***************
          end
          
          %*********** Here evaluate goodness of fit of VonMises
          fit1 = Info.SacOn.OriTune{1}.Afit;  % constrained same pref
          afit1 = Info.SacOn.OriTune{1};
          oval1 = Info.SacOn.OriTune{1}.OriVals;
          %***
          fit2 = Info.SacOn.OriTune{2}.Afit;
          afit2 = Info.SacOn.OriTune{2};
          oval2 = Info.SacOn.OriTune{2}.OriVals;
          %***
          sfit1 = Info.StimOn.OriTune{1}.Afit;
          sfit2 = Info.StimOn.OriTune{2}.Afit;
          %****************
          %*** could do based on R2 of the fit instead, using
          %*** raw ori tune values against the vonmises fits
          %*********
          % throw out if kappa parameter at its limit (bad fit?)
          if (1)
              %*** code for computations
              acur = VonMises.vonmises(oval1*(pi/180),fit1.mu)';
              bcur = VonMises.vonmises(oval2*(pi/180),fit2.mu)';
              %******
              Ea = nansum( (afit1.Mu - acur) .^ 2);
              Eta = nansum( (afit1.Mu - nanmean(afit1.Mu)) .^ 2);
              Eb = nansum( (afit2.Mu - bcur) .^ 2);
              Etb = nansum( (afit2.Mu - nanmean(afit2.Mu)) .^ 2);
              %*****
              R2a = (Eta - Ea)/Eta;
              R2b = (Etb - Eb)/Etb;
              R2 = ((Eta + Etb)-(Ea+Eb))/(Eta+Etb);  % total fit, allows for
                                                     % null tuning in one
                                                     % cond if strong other
              
                                                     
              aimic = ((Info.AUC_SacOn.MIC{2,1}(1)-Info.AUC_SacOn.MIC{2,1}(2))/...
                   (Info.AUC_SacOn.MIC{2,1}(1)+Info.AUC_SacOn.MIC{2,1}(2)));
              %if ( (aimic > 0.08) & (R2 > 0.7) )
              if (0)    
                  fname
                  
              %******* Make a nice single unit plot (Shanna thesis)
              % if (0) %(strcmp(fname,'MT_110719_U2'))
                  
                  figure(10); hold off;
                
                  vmax = max( [max(afit1.Atune.mu+(2*afit1.Atune.sem)),max(afit2.Atune.mu+(2*afit2.Atune.sem)),...
                               max(afit1.Mu+2*afit1.Sem), max(afit2.Mu+2*afit2.Sem)] );
                  vmax = vmax * 1.1;  
                  aa = [0 340 340 0 0];
                  bb = [vmax vmax 0 0 vmax];
                  fill(aa,bb,[1,1,1],'FaceAlpha',1.0,'Linestyle','none'); hold on;
                  
                  sumtune = afit1.Atune.mu + afit2.Atune.mu;
                  zp = find( sumtune == max(sumtune));
                  %*******
                  axx = afit1.OriVals';
                  aN = length(axx);
                  amid = floor(aN/2);
                  ashift = zp(1)-amid;
                  if (ashift == 0)
                      zaxx = [1:aN];
                  else
                    if (ashift < 0)
                      zaxx = [(aN+ashift):aN,1:(aN+ashift-1)]; % -ashift:aN,1:(-ashift-1)];
                    else
                      zaxx = [ashift:aN,1:(ashift-1)];    
                    end
                  end
                  %*******
                  disp('Plotting ..');
                  disp(sprintf('Iso %d Shift %d',iso,ashift));
                  
                  plot(oval1,afit1.Mu(zaxx),'r.','Markersize',10); hold on;
                  for ii = 1:length(oval1)
                      plot([oval1(ii),oval1(ii)],[afit1.Mu(zaxx(ii))+2*afit1.Sem(zaxx(ii)),...
                                         afit1.Mu(zaxx(ii))-2*afit1.Sem(zaxx(ii))],'r-');
                  end
                  %******
                  plot(oval2,afit2.Mu(zaxx),'b.','Markersize',10); hold on;
                  for ii = 1:length(oval2)
                      plot([oval2(ii),oval2(ii)],[afit2.Mu(zaxx(ii))+2*afit2.Sem(zaxx(ii)),...
                                         afit2.Mu(zaxx(ii))-2*afit2.Sem(zaxx(ii))],'b-');
                  end
                  %********
                  
                  colo = [1,0,0];
                  % tta = [afit1.OriVals' fliplr(afit1.OriVals')];
                  % uua = [(afit1.Atune.mu+(2*afit1.Atune.sem)) ...
                  %                fliplr((afit1.Atune.mu-(2*afit1.Atune.sem)))];
                  % fill(tta,uua,colo,'FaceAlpha',0.2,'Linestyle','none'); hold on;
                  plot(afit1.OriVals',afit1.Atune.mu(zaxx),'k-','Color',colo,'Linewidth',2);
                  plot(afit1.OriVals',afit1.Atune.mu(zaxx)+(2*afit2.Atune.sem(zaxx)),'k-','Color',colo);
                  plot(afit1.OriVals',afit1.Atune.mu(zaxx)-(2*afit2.Atune.sem(zaxx)),'k-','Color',colo);
  
                  %******
                  colo = [0,0,1];
                  %tta = [afit2.OriVals' fliplr(afit2.OriVals')];
                  %uua = [(afit2.Atune.mu+(2*afit2.Atune.sem)) ...
                  %                fliplr((afit2.Atune.mu-(2*afit2.Atune.sem)))];
                  %fill(tta,uua,colo,'FaceAlpha',0.2,'Linestyle','none'); hold on;
                  plot(afit2.OriVals',afit2.Atune.mu(zaxx),'k-','Color',colo,'Linewidth',2);
                  plot(afit2.OriVals',afit2.Atune.mu(zaxx)+(2*afit2.Atune.sem(zaxx)),'k-','Color',colo);
                  plot(afit2.OriVals',afit2.Atune.mu(zaxx)-(2*afit2.Atune.sem(zaxx)),'k-','Color',colo);
                  %******
  
                  axis tight;
                  V = axis;
                  axis([-30 340 0 vmax]);
                  axis off;
                  plot([180,340],[0,0],'k-','Linewidth',2);
                  text(250,-(vmax*0.04),'180 degs','Fontsize',18);
                  vscal = floor(floor(vmax*0.5)/10)*10;
                  plot([-20,-20],[0,vscal],'k-','Linewidth',2);
                  text(-40,vscal*0.3,sprintf('%2d sp/s',vscal),'Rotation',90,'Fontsize',18);
                  
                 % title(sprintf('%s R2:%5.3f Base,Gain,Width(%+4.2f,%+4.2f,%+4.2f)Iso(%d)',fname,R2,...
                 %        (AttBase-IgnBase)/(AttBase+IgnBase),...
                 %        (AttGain-IgnGain)/(AttGain+IgnGain),...
                 %        (AttWid-IgnWid)/(AttWid+IgnWid),iso));
                 input('check');
              end

          end
     
          
          %**** Goodness of fit criterion for curve analyses
          R2Pass = 1;
          if (1)
            if (R2 < Rthresh)
              R2Pass = 0;
            end
          end
          R2List = [R2List ; R2Pass];  % include for curve analysis or not
     
          if (strcmp(fname,ExampleName))
              ExampleInfo = Info;
              ExampleInd = length(R2List);  %record its index
          end  
          
          %*********** Here evaluate goodness of fit of VonMises
          zfit1 = Info.SacOn.FeaTune{1}.Afit;
          zafit1 = Info.SacOn.FeaTune{1};
          zoval1 = Info.SacOn.FeaTune{1}.OriVals;
          zfit2 = Info.SacOn.FeaTune{2}.Afit;
          zafit2 = Info.SacOn.FeaTune{2};
          zoval2 = Info.SacOn.FeaTune{2}.OriVals;        
          
          %****** non-parametric smoothing of tuning curves
          %******* compute rank ordered tuning curve
          LO = length(afit1.Mu);
          NPSMOOTH = 1;
          aa = circsmooth(afit1.Mu,NPSMOOTH);
          uu = circsmooth(afit2.Mu,NPSMOOTH);
          sa = circsmooth(afit1.Sem,NPSMOOTH);
          su = circsmooth(afit2.Sem,NPSMOOTH);
          tot = aa + uu;
          theta = 0:(LO-1);
          theta = theta*(2*pi)/LO;
          wvec = sum( tot .* exp(j*theta) );
          wvec = wvec/abs(wvec);
          awvec = sum( aa .* exp(j*theta) );
          awvec = awvec/abs(awvec);
          uwvec = sum( uu .* exp(j*theta) );
          uwvec = uwvec/abs(uwvec);
          % force equal weighting of center between conditions
          if (centerwid) 
              wvec = 0.5*(awvec+uwvec);
              wvec = wvec / abs(wvec);
          end
          %*****
          rval = [];
          for ii = 1:length(oval1)
                    ang = ((ii-1)*2*pi/length(oval1));
                    dot = (real(wvec)*cos(ang) + imag(wvec)*sin(ang));
                    rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
          end
          yrvals = sortrows(rval,2);
          rord = yrvals(:,1);
          %*******
          rval = [];
          for ii = 1:length(oval1)
                    ang = ((ii-1)*2*pi/length(oval1));
                    dot = (real(awvec)*cos(ang) + imag(awvec)*sin(ang));
                    rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
          end
          yrvals = sortrows(rval,2);
          arord = yrvals(:,1);
          %*******
          rval = [];
          for ii = 1:length(oval1)
                    ang = ((ii-1)*2*pi/length(oval1));
                    dot = (real(uwvec)*cos(ang) + imag(uwvec)*sin(ang));
                    rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
          end
          yrvals = sortrows(rval,2);
          urord = yrvals(:,1);
          %*******
          ratt = aa(rord);  %smoothnp(afit1.Mu(rord),NPSMOOTH);
          sratt = sa(rord); %smoothnp(afit1.Sem(rord),NPSMOOTH);
          rutt = uu(rord);  %smoothnp(afit2.Mu(rord),NPSMOOTH);
          srutt = su(rord); %smoothnp(afit2.Sem(rord),NPSMOOTH);
          [ampw,umpw,abasew,ubasew,aw,uw] = compute_NP_params(ratt,rutt);
          zratt = smoothnp(zafit1.Mu(rord),NPSMOOTH);
          zsratt = smoothnp(zafit1.Sem(rord),NPSMOOTH);
          zrutt = smoothnp(zafit2.Mu(rord),NPSMOOTH);
          zsrutt = smoothnp(zafit2.Sem(rord),NPSMOOTH);
          
          %****** non-parametric smoothing of tuning curves
          %******* compute rank ordered tuning curve
          rafit1 = Info.SacOn.RandOriTune{1};
          rafit2 = Info.SacOn.RandOriTune{2};
          LO = length(rafit1.Mu);
          NPSMOOTH = 1;
          aa = circsmooth(rafit1.Mu,NPSMOOTH);
          uu = circsmooth(rafit2.Mu,NPSMOOTH);
          tot = aa + uu;
          theta = 0:(LO-1);
          theta = theta*(2*pi)/LO;
          wvec = sum( tot .* exp(j*theta) );
          wvec = wvec/abs(wvec);
          awvec = sum( aa .* exp(j*theta) );
          awvec = awvec/abs(awvec);
          uwvec = sum( uu .* exp(j*theta) );
          uwvec = uwvec/abs(uwvec);
          % force equal weighting of center between conditions
          if (centerwid) 
              wvec = 0.5*(awvec+uwvec);
              wvec = wvec / abs(wvec);
          end
          %*****
          rval = [];
          for ii = 1:length(oval1)
                    ang = ((ii-1)*2*pi/length(oval1));
                    dot = (real(wvec)*cos(ang) + imag(wvec)*sin(ang));
                    rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
          end
          yrvals = sortrows(rval,2);
          rord = yrvals(:,1);
          %*******
          rratt = aa(rord);  %smoothnp(afit1.Mu(rord),NPSMOOTH);
          rrutt = uu(rord);  %smoothnp(afit2.Mu(rord),NPSMOOTH);
          [rampw,rumpw,rabasew,rubasew,raw,ruw] = compute_NP_params(rratt,rrutt);
          
          %******* compute width with error bars here ***
          aamps = []; uamps = [];
          abases = []; ubases = [];
          awids = []; uwids = [];
          for jk = 1:size(Info.SacOn.OriTune{1}.JACKTUNE,1)
              jackaa = Info.SacOn.OriTune{1}.JACKTUNE(jk,:); 
              jackuu = Info.SacOn.OriTune{2}.JACKTUNE(jk,:);
              %******
              LO = length(jackaa);
              NPSMOOTH = 1;
              aa = circsmooth(jackaa,NPSMOOTH);
              uu = circsmooth(jackuu,NPSMOOTH);
              tot = aa + uu;
              theta = 0:(LO-1);
              theta = theta*(2*pi)/LO;
              wvec = sum( tot .* exp(j*theta) );
              wvec = wvec/abs(wvec);
              %***************
              if (centerwid) 
                  awvec = sum( aa .* exp(j*theta) );
                  awvec = awvec/abs(awvec);
                  uwvec = sum( uu .* exp(j*theta) );
                  uwvec = uwvec/abs(uwvec);
                  %*******
                  wvec = 0.5*(awvec+uwvec);
                  wvec = wvec / abs(wvec);
              end
              %*********
              Wvec2 = wvec;
              ango = angle(wvec);
              if (ango < 0)
                 ango = ango + 2*pi;
              end
              %*****
              rval = [];
              for ii = 1:length(oval1)
                        ang = ((ii-1)*2*pi/length(oval1));
                        dot = (real(Wvec2)*cos(ang) + imag(Wvec2)*sin(ang));
                        rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
              end
              yrvals = sortrows(rval,2);
              rord = yrvals(:,1);
              rratt = jackaa(rord);  %smoothnp(afit1.Mu(rord),NPSMOOTH);
              rrutt = jackuu(rord);  %smoothnp(afit2.Mu(rord),NPSMOOTH);
              [ampa,ampu,basea,baseu,awid,uwid] = compute_NP_params(rratt,rrutt);
              %*****
              aamps = [aamps ; ampa];
              uamps = [uamps ; ampu];
              abases = [abases ; basea];
              ubases = [ubases ; baseu];
              awids = [awids ; awid];
              uwids = [uwids ; uwid];
          end
          %***** Jacknife estimates on widths (which units are significant)
          NW = length(awids)*(1+2*NPSMOOTH);
          ampsw = nanstd(aamps)*sqrt((NW-1));
          umpsw = nanstd(uamps)*sqrt((NW-1));
          abasesw = nanstd(abases)*sqrt((NW-1));
          ubasesw = nanstd(ubases)*sqrt((NW-1));
          asw = nanstd(awids)*sqrt((NW-1));
          usw = nanstd(uwids)*sqrt((NW-1));
          %***************
 
          %******* compute width with error bars here ***
          raamps = []; ruamps = [];
          rabases = []; rubases = [];
          rawids = []; ruwids = [];
          for jk = 1:size(Info.SacOn.RandOriTune{1}.JACKTUNE,1)
              jackaa = Info.SacOn.RandOriTune{1}.JACKTUNE(jk,:); 
              jackuu = Info.SacOn.RandOriTune{2}.JACKTUNE(jk,:);
              %******
              LO = length(jackaa);
              NPSMOOTH = 1;
              aa = circsmooth(jackaa,NPSMOOTH);
              uu = circsmooth(jackuu,NPSMOOTH);
              tot = aa + uu;
              theta = 0:(LO-1);
              theta = theta*(2*pi)/LO;
              wvec = sum( tot .* exp(j*theta) );
              wvec = wvec/abs(wvec);
              if (centerwid) 
                  awvec = sum( jackaa .* exp(j*theta) );
                  awvec = awvec/abs(awvec);
                  uwvec = sum( jackuu .* exp(j*theta) );
                  uwvec = uwvec/abs(uwvec);
                  %*******
                  wvec = 0.5*(awvec+uwvec);
                  wvec = wvec / abs(wvec);
              end
              Wvec2 = wvec;
              ango = angle(wvec);
              if (ango < 0)
                 ango = ango + 2*pi;
              end
              %*****
              rval = [];
              for ii = 1:length(oval1)
                        ang = ((ii-1)*2*pi/length(oval1));
                        dot = (real(Wvec2)*cos(ang) + imag(Wvec2)*sin(ang));
                        rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
              end
              yrvals = sortrows(rval,2);
              rord = yrvals(:,1);
              rratt = jackaa(rord);  %smoothnp(afit1.Mu(rord),NPSMOOTH);
              rrutt = jackuu(rord);  %smoothnp(afit2.Mu(rord),NPSMOOTH);
              [rampa,rampu,rbasea,rbaseu,rawid,ruwid] = compute_NP_params(rratt,rrutt);
              %*****
              raamps = [raamps ; rampa];
              ruamps = [ruamps ; rampu];
              rabases = [rabases ; rbasea];
              rubases = [rubases ; rbaseu];
              rawids = [rawids ; rawid];
              ruwids = [ruwids ; ruwid];
          end
          %***** Jacknife estimates on widths (which units are significant)
          NW = length(rawids)*(1+2*NPSMOOTH);
          rampsw = nanstd(raamps)*sqrt((NW-1));
          rumpsw = nanstd(ruamps)*sqrt((NW-1));
          rabasesw = nanstd(rabases)*sqrt((NW-1));
          rubasesw = nanstd(rubases)*sqrt((NW-1));
          rasw = nanstd(rawids)*sqrt((NW-1));
          rusw = nanstd(ruwids)*sqrt((NW-1));
          %***************
          
          %****************
          %*** could do based on R2 of the fit instead, using
          %*** raw ori tune values against the vonmises fits
          %*********
          % throw out if kappa parameter at its limit (bad fit?)
          if (1)
              %*** code for computations
              acur = VonMises.vonmises(oval1*(pi/180),zfit1.mu)';
              bcur = VonMises.vonmises(oval2*(pi/180),zfit2.mu)';
              %******
              Ea = nansum( (zafit1.Mu - acur) .^ 2);
              Eta = nansum( (zafit1.Mu - nanmean(zafit1.Mu)) .^ 2);
              Eb = nansum( (zafit2.Mu - bcur) .^ 2);
              Etb = nansum( (zafit2.Mu - nanmean(zafit2.Mu)) .^ 2);
              %*****
              zR2a = (Eta - Ea)/Eta;
              zR2b = (Etb - Eb)/Etb;
              zR2 = ((Eta + Etb)-(Ea+Eb))/(Eta+Etb);  % total fit, allows for
                                                     % null tuning in one
                                                     % cond if strong other
              %*********
              if (0) % (k > N1) & (R2 > 0.5)
                  h10 = figure(10); 
                  set(h10,'Position',[100 100 800 400]);
                  subplot('Position',[0.1 0.15 0.35 0.70]); hold off;
                  plot(oval1,zafit1.Mu,'ro-'); hold on;
                  for ii = 1:length(oval1)
                      plot([oval1(ii),oval1(ii)],[zafit1.Mu(ii)+2*zafit1.Sem(ii),...
                                         zafit1.Mu(ii)-2*zafit1.Sem(ii)],'r-');
                  end
                  plot(oval1,acur,'r--','LineWidth',2);
                  plot(oval2,zafit2.Mu,'bo-'); hold on;
                  for ii = 1:length(oval2)
                      plot([oval2(ii),oval2(ii)],[zafit2.Mu(ii)+2*zafit2.Sem(ii),...
                                         zafit2.Mu(ii)-2*zafit2.Sem(ii)],'b-');
                  end
                  plot(oval2,bcur,'b--','LineWidth',2);
                  title(sprintf('R2:(%5.3f,%5.3f)Tot(%5.3f)',zR2a,zR2b,zR2));
                  %*****
                  subplot('Position',[0.6 0.15 0.35 0.70]); hold off;
                  plot(1:length(oval1),ratt,'r.-'); hold on;
                  for ii = 1:length(oval1)
                      plot([ii,ii],[ratt(ii)+2*sratt(ii),...
                                    ratt(ii)-2*sratt(ii)],'r-');
                  end
                  plot(1:length(oval1),rutt,'b.-'); hold on;
                  for ii = 1:length(oval1)
                      plot([ii,ii],[rutt(ii)+2*srutt(ii),...
                                    rutt(ii)-2*srutt(ii)],'b-');
                  end
                  plot(1:length(oval1),rrank,'k.-'); hold on;
                  for ii = 1:length(oval1)
                      plot([ii,ii],[rrank(ii)+2*srank(ii),...
                                    rrank(ii)-2*srank(ii)],'k-');
                  end
                  plot(1:length(oval1),zratt,'ro-'); hold on;
                  for ii = 1:length(oval1)
                      plot([ii,ii],[zratt(ii)+2*zsratt(ii),...
                                    zratt(ii)-2*zsratt(ii)],'r-');
                  end
                  plot(1:length(oval1),zrutt,'bo-'); hold on;
                  for ii = 1:length(oval1)
                      plot([ii,ii],[zrutt(ii)+2*zsrutt(ii),...
                                    zrutt(ii)-2*zsrutt(ii)],'b-');
                  end
                  xlabel('Rank Order');
                  ylabel('Firing Rate');
                  input('check');
              end
          end
          
          %**** Goodness of fit criterion for curve analyses
          zR2Pass = 1;
          if (1)
            if (zR2 < Rthresh)
              zR2Pass = 0;
            end
          end
          zR2List = [zR2List ; zR2Pass];  % include for curve analysis or not
           
          %*************
          OriVals = Info.OriTune.OriVals;
          LockWin = Info.SacOn.LockWin;
          LockInt = Info.SacOn.LockInt;
          Event = Info.SacOn.Event;
          LockWin2 = Info.StimOn.LockWin;
          LockInt2 = Info.StimOn.LockInt;
          Event2 = Info.StimOn.Event;
          %******
          NDepth = [NDepth; Info.depth];
          Layer = [Layer; Info.layer];
          CSD = [CSD; Info.CSD];
          Duration = [Duration ; Info.duration];
          %*******
          DSI{1} = [DSI{1}; [Info.OriTune.DSI]]; %%New AB
          VONFIT{1} = [VONFIT{1} ; [animal, iso, fit1.mu, fit2.mu]];      
          VONFITSEM{1} = [VONFITSEM{1} ; [animal, iso, fit1.sem, fit2.sem]];      
          VONFIT{2} = [VONFIT{2} ; [animal, iso, sfit1.mu, sfit2.mu]];      
          VONFITSEM{2} = [VONFITSEM{2} ; [animal, iso, sfit1.sem, sfit2.sem]];      
          ZVONFIT{1} = [ZVONFIT{1} ; [animal, iso, zfit1.mu, zfit2.mu]];      
          ZVONFITSEM{1} = [ZVONFITSEM{1} ; [animal, iso, zfit1.sem, zfit2.sem]];      

          %********
          NPFIT{1} = [NPFIT{1} ; ratt];  % non-parametric curve fits
          NPFIT{2} = [NPFIT{2} ; rutt];
          NPFEAT{1} = [NPFEAT{1} ; zratt];  % non-parametric curve fits (feature)
          NPFEAT{2} = [NPFEAT{2} ; zrutt];
          %*******
          NPAMP = [NPAMP ; [animal, ampw, ampsw, umpw, umpsw]];
          NPBASE = [NPBASE ; [animal, abasew, abasesw, ubasew, ubasesw]];
          NPWID = [NPWID ; [animal, aw, asw, uw, usw]];  % non-param wid with sems
          NPMIN = [NPMIN ; minrate];  % min firing rate of two conditions
          
          RNPAMP = [RNPAMP ; [animal, rampw, rampsw, rumpw, rumpsw]];
          RNPBASE = [RNPBASE ; [animal, rabasew, rabasesw, rubasew, rubasesw]];
          RNPWID = [RNPWID ; [animal, raw, rasw, ruw, rusw]];  % non-param wid with sems
          
          %****** measure LOCATION BIAS TOO
          Na = length(Info.AttList);   % Ephys trials with target in RF
          Nb = length(Info.IgnAList) + length(Info.IgnBList); 
          LOCBIAS = [LOCBIAS ; ((Na-Nb)/(Na+Nb))];
          
          %***** Baseline fits of Von Mises
          AttBase = fit1.mu(1);
          AttBaseSem = fit1.sem(1);
          IgnBase = fit2.mu(1);
          IgnBaseSem = fit2.sem(1);
          %***** Baseline fits of Von Mises
          AttBase2 = sfit1.mu(1);
          AttBaseSem2 = sfit1.sem(1);
          IgnBase2 = sfit2.mu(1);
          IgnBaseSem2 = sfit2.sem(1);
          %***** Gain fits of Von Mises
          AttGain = fit1.mu(2);
          AttGainSem = fit1.sem(2);
          IgnGain = fit2.mu(2);
          IgnGainSem = fit2.sem(2);
          %***** Gain fits of Von Mises
          AttGain2 = sfit1.mu(2);
          AttGainSem2 = sfit1.sem(2);
          IgnGain2 = sfit2.mu(2);
          IgnGainSem2 = sfit2.sem(2);
          % estimate the actual half-widths, not kappa
          %******* compute half-width
          [AttWid,IgnWid] = compute_half_width(fit1.mu,fit2.mu);
          [AttWid2,IgnWid2] = compute_half_width(sfit1.mu,sfit2.mu);
          %**** use sampling to estimate half-width error bars
          AW = []; IW = [];
          AW2 = []; IW2 = [];
          for it = 1:100
               mu1 = fit1.mu + ([0,0,randn,0] .* fit1.sem);
               mu2 = fit2.mu + ([0,0,randn,0] .* fit2.sem);
               [aw,iw] = compute_half_width(mu1,mu2);
               AW = [AW ; aw];
               IW = [IW ; iw];
               %*****
               mu1 = sfit1.mu + ([0,0,randn,0] .* sfit1.sem);
               mu2 = sfit2.mu + ([0,0,randn,0] .* sfit2.sem);
               [aw2,iw2] = compute_half_width(sfit1.mu,sfit2.mu);
               AW2 = [AW2 ; aw2];
               IW2 = [IW2 ; iw2];
               %****
          end
          AttWidSem = nanstd(AW);
          IgnWidSem = nanstd(IW);
          AttWidSem2 = nanstd(AW2);
          IgnWidSem2 = nanstd(IW2);
          %******** min is 0 and max is 180  
          
          %***** Baseline fits of Von Mises
          zAttBase = zfit1.mu(1);
          zAttBaseSem = zfit1.sem(1);
          zIgnBase = zfit2.mu(1);
          zIgnBaseSem = zfit2.sem(1);
          %***** Gain fits of Von Mises
          zAttGain = zfit1.mu(2);
          zAttGainSem = zfit1.sem(2);
          zIgnGain = zfit2.mu(2);
          zIgnGainSem = zfit2.sem(2);
          % estimate the actual half-widths, not kappa
          %******* compute half-width
          [zAttWid,zIgnWid] = compute_half_width(zfit1.mu,zfit2.mu);
          zAttWidSem = NaN;
          zIgnWidSem = NaN;
          %******** min is 0 and max is 180  
          
          
          %***** Angular Std fits of Von Mises
          fit1 = Info.SacOn.OriTune{1};
          fit2 = Info.SacOn.OriTune{2};
          sfit1 = Info.SacOn.OriTune{1};
          sfit2 = Info.SacOn.OriTune{2};
          %******
          AttAng = fit1.DSI; %AngStd;
          AttAngSem = fit1.DSI_SEM; %ANGSTD_SEM;
          IgnAng = fit2.DSI; %AngStd; 
          IgnAngSem = fit2.DSI_SEM; %ANGSTD_SEM; 
          %***** Angular Std fits of Von Mises
          AttAng2 = sfit1.DSI; %AngStd;
          AttAngSem2 = sfit1.DSI_SEM; %ANGSTD_SEM;
          IgnAng2 = sfit2.DSI; %AngStd; 
          IgnAngSem2 = sfit2.DSI_SEM; %ANGSTD_SEM;
          
          %*****
          %****** Store AUC in main interval (1 is -50 to 30)
          klock = 1;  % -50 to 30 ms
          nsmo = 2;  % bin smoothing of 1 in circular mutual info if this is 2
          nsmo2 = 1;
          ATTAUC = [ATTAUC ; Info.AUC_SacOn.AUC{klock}(1)];
          IGNAUC = [IGNAUC ; Info.AUC_SacOn.AUC{klock}(2)];
          %****
          SacOn_TT = Info.SacOn.LockTT{1};
          ztt = find( (SacOn_TT >= LockInt(1)) & (SacOn_TT < LockInt(2)));
          SacOn_Arate = nanmean(Info.SacOn.LockUU{1}(ztt));
          SacOn_Irate = nanmean(Info.SacOn.LockUU{2}(ztt));
          %*****
          StimOn_TT = Info.StimOn.LockTT{1};
          ztt = find( (StimOn_TT >= LockInt2(1)) & (StimOn_TT < LockInt2(2)));
          StimOn_Arate = nanmean(Info.StimOn.LockUU{1}(ztt));
          StimOn_Irate = nanmean(Info.StimOn.LockUU{2}(ztt));
          %**** 
          FR{1} = [FR{1} ; [animal, iso, SacOn_Arate, 0, SacOn_Irate,0,tungsten]];
          AUC{1} = [AUC{1} ; [animal, iso, Info.AUC_SacOn.AUC{klock}(1), 0, Info.AUC_SacOn.AUC{klock}(2),0,tungsten]];
          MI{1} = [MI{1} ; [animal, iso, Info.AUC_SacOn.MI{klock}(1), 0, Info.AUC_SacOn.MI{klock}(2), 0, tungsten]];
          MIC{1} = [MIC{1} ; [animal, iso, Info.AUC_SacOn.MIC{nsmo,klock}(1), 0, Info.AUC_SacOn.MIC{nsmo,klock}(2), 0, tungsten]];
          MIC2{1} = [MIC2{1} ; [animal, iso, Info.AUC_SacOn.MIC{nsmo2,klock}(1), 0, Info.AUC_SacOn.MIC{nsmo2,klock}(2), 0, tungsten]];
          % FF{1} = [FF{1} ; [animal, iso, Info.AUC_SacOn.Fano{klock}(1), 0, Info.AUC_SacOn.Fano{klock}(2), 0, tungsten]];
          FF{1} = [FF{1} ; [animal, iso, Info.AUC_SacOn.RFano{klock}(1), 0, Info.AUC_SacOn.RFano{klock}(2), 0, tungsten]];
          FS{1} = [FS{1} ; [animal, iso, Info.AUC_SacOn.FanoSlope{klock}(1), 0, Info.AUC_SacOn.FanoSlope{klock}(2), 0, tungsten]];
          %****
          FR{2} = [FR{2} ; [animal, iso, StimOn_Arate, 0, StimOn_Irate,0,tungsten]];
          AUC{2} = [AUC{2} ; [animal, iso, Info.AUC_StimOn.AUC{klock}(1), 0, Info.AUC_StimOn.AUC{klock}(2),0,tungsten]];
          MI{2} = [MI{2} ; [animal, iso, Info.AUC_StimOn.MI{klock}(1), 0, Info.AUC_StimOn.MI{klock}(2), 0, tungsten]];
          MIC{2} = [MIC{2} ; [animal, iso, Info.AUC_StimOn.MIC{nsmo,klock}(1), 0, Info.AUC_StimOn.MIC{nsmo,klock}(2), 0, tungsten]];
          MIC2{2} = [MIC2{2} ; [animal, iso, Info.AUC_StimOn.MIC{nsmo2,klock}(1), 0, Info.AUC_StimOn.MIC{nsmo2,klock}(2), 0, tungsten]];
          %FF{2} = [FF{2} ; [animal, iso, Info.AUC_StimOn.Fano{klock}(1), 0, Info.AUC_StimOn.Fano{klock}(2), 0, tungsten]];
          FF{2} = [FF{2} ; [animal, iso, Info.AUC_StimOn.RFano{klock}(1), 0, Info.AUC_StimOn.RFano{klock}(2), 0, tungsten]];
          FS{2} = [FS{2} ; [animal, iso, Info.AUC_StimOn.FanoSlope{klock}(1), 0, Info.AUC_StimOn.FanoSlope{klock}(2), 0, tungsten]];
          %******
          VONBASE{1} = [VONBASE{1} ; [animal, iso, AttBase, AttBaseSem, IgnBase, IgnBaseSem,tungsten]];
          VONBASE{2} = [VONBASE{2} ; [animal, iso, AttBase2, AttBaseSem2, IgnBase2, IgnBaseSem2,tungsten]];
          VONGAIN{1} = [VONGAIN{1} ; [animal, iso, AttGain, AttGainSem, IgnGain, IgnGainSem,tungsten]];
          VONGAIN{2} = [VONGAIN{2} ; [animal, iso, AttGain2, AttGainSem2, IgnGain2, IgnGainSem2,tungsten]];
          VONKAP{1} = [VONKAP{1} ; [animal, iso, AttWid, AttWidSem, IgnWid, IgnWidSem,tungsten]];
          VONKAP{2} = [VONKAP{2} ; [animal, iso, AttWid2, AttWidSem2, IgnWid2, IgnWidSem2,tungsten]];
          %********
          ZVONBASE{1} = [ZVONBASE{1} ; [animal, iso, zAttBase, zAttBaseSem, zIgnBase, zIgnBaseSem,tungsten]];
          ZVONGAIN{1} = [ZVONGAIN{1} ; [animal, iso, zAttGain, zAttGainSem, zIgnGain, zIgnGainSem,tungsten]];
          ZVONKAP{1} = [ZVONKAP{1} ; [animal, iso, zAttWid, zAttWidSem, zIgnWid, zIgnWidSem,tungsten]];
          %********
          RATE_tt{1} = [RATE_tt{1} ; Info.SacOn.LockTT{1} ];
          RATE_att{1} = [RATE_att{1} ; Info.SacOn.LockUU{1} ];
          RATE_ign{1} = [RATE_ign{1} ; Info.SacOn.LockUU{2} ];
          RATE_tt{2} = [RATE_tt{2} ; Info.StimOn.LockTT{1} ];
          RATE_att{2} = [RATE_att{2} ; Info.StimOn.LockUU{1} ];
          RATE_ign{2} = [RATE_ign{2} ; Info.StimOn.LockUU{2} ];
          RATE_Anim = [RATE_Anim ; animal];
          %*******
          if (0)  % Sanity check the fits
            if ((AttWid2/IgnWid2) > 2)
              figure(10); hold off;
              [AttBase2,IgnBase2]
              [AttGain2,IgnGain2]
              [AttWid2,IgnWid2]
              Afit = Info.SacWin2.OriTune{1}.Afit
              Ufit = Info.SacWin2.OriTune{2}.Afit
              x = 0:0.01:(2*pi); 
              aa = VonMises.vonmises(x,Afit.mu);
              uu = VonMises.vonmises(x,Ufit.mu);
              plot(x,aa,'r-'); hold on;
              plot(x,uu,'b-');
              input('check');
            end
            clear Info;
          end
      end
      %****
  end
  %********* finished processing, then store environment to reload next
  %********* for the next time ***********
  save(Environment_Preload)
  %*********
end

%****** report on cells included for each animal
disp(sprintf('MTC:  total(%d)(%5.2f vis include)',Etot,100*((Eiso+Emu)/Etot)));
disp(sprintf('MTC:  isolated (%d)  multi-unit(%d)',Eiso,Emu));
disp(sprintf('MTC:  include total(%d) tuned(%d)(%5.2f dsi sig) Incrate(%d)',...
             (Eiso+Emu),Etuned,100*(Etuned/(Eiso+Emu)),Eincrate));
disp(sprintf('MTC:  final included (single %d, multi %d',Esingle,Emulti));         
disp('********');
disp(sprintf('MT:  total(%d)(%5.2f vis include)',Mtot,100*((Miso+Mmu)/Mtot)));
disp(sprintf('MT:  isolated (%d)  multi-unit(%d)',Miso,Mmu));
disp(sprintf('MT:  include total(%d) tuned(%d)(%5.2f dsi sig) Incrate(%d)',...
             (Miso+Mmu),Mtuned,100*(Mtuned/(Miso+Mmu)),Mincrate));
disp(sprintf('MT:  Mu tungsten (%d)',Mtungsten));
disp(sprintf('MT:  final included (single %d, multi %d',Msingle,Mmulti));         
%**************************

%*********  make nice figure for paper/thesis
if (1)  
  Mcolo = [0,0.6,0.6];
  Ecolo = [0.3,0.7,0.3];
  hf = figure;
  set(hf,'Position',[100 100 1200 800]);
  %***** first, plot the PSTH per animal on the left
  for zkk = 1:2

    %*****
    if (zkk == 1)
      zz1 = find( RATE_Anim == 1);  % MTC 1
      % zz1 = narrow;
      subplot('Position',[0.10 0.60 0.4 0.35]);
    else
      zz1 = find( RATE_Anim == 2);  % MTC 1
      % zz1 = broad;
      subplot('Position',[0.10 0.10 0.4 0.35]);
    end
    tt1 = nanmean(RATE_tt{1}(zz1,:));
    tt2 = nanmean(RATE_tt{2}(zz1,:));
    att = cell(2,2);
    ign = cell(2,2);
    for kk = 1:2
        arate = [];
        irate = [];
        maxos = [];
        for zzk = 1:size(zz1,1)
          ma = [];
          for uk = 1:2  % do over both conditions so we can keep on same scale
           ma = [ma max(RATE_att{uk}(zz1(zzk),:))];
           ma = [ma max(RATE_ign{kk}(zz1(zzk),:))];
          end
          maxo = max(ma);
          maxos = [maxos maxo];
          arate = [arate ; (RATE_att{kk}(zz1(zzk),:)/maxo)];
          irate = [irate ; (RATE_ign{kk}(zz1(zzk),:)/maxo)];
        end
        scal = 1; % mean(maxos);
        att{kk,1} = scal * nanmean(arate);
        att{kk,2} = scal * nanstd(arate)/sqrt(size(arate,1));
        utt{kk,1} = scal * nanmean(irate);
        utt{kk,2} = scal * nanstd(irate)/sqrt(size(irate,1));
    end
    %*****
    if (zkk == 1) % set max based from narrow max rates
       % vmax = max(att{2,1}+(2*att{2,2}))*1.1;
       vmax = 0.85;
    end
    
    aa = [-0.05 0.25 0.25 -0.05 -0.05];
    bb = [vmax vmax 0 0 vmax];
    fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
    kk = 1;
    zz = find( (tt1 >= -0.05) & (tt1 < 0.05) );
    bd = 0.20;
    aa = [(bd+tt1(zz)) fliplr(bd+tt1(zz))];
    bb = [(att{kk,1}(zz)+2*att{kk,2}(zz)) fliplr(att{kk,1}(zz)-2*att{kk,2}(zz))];
    fill(aa,bb,[1,0,0],'Linestyle','none','FaceAlpha',0.3); hold on;
    plot(bd+tt1(zz),att{kk,1}(zz),'r-','Linewidth',2); hold on;
    bb = [(utt{kk,1}(zz)+2*utt{kk,2}(zz)) fliplr(utt{kk,1}(zz)-2*utt{kk,2}(zz))];
    fill(aa,bb,[0,0,1],'Linestyle','none','FaceAlpha',0.3); hold on;
    plot(bd+tt1(zz),utt{kk,1}(zz),'b-','Linewidth',2); hold on;
    %******
    kk = 2;
    zz = find( (tt2 >= -0.025) & (tt2 < 0.10) );
    bd2 = 0.0;
    aa = [(bd2+tt2(zz)) fliplr(bd2+tt2(zz))];
    bb = [(att{kk,1}(zz)+2*att{kk,2}(zz)) fliplr(att{kk,1}(zz)-2*att{kk,2}(zz))];
    fill(aa,bb,[1,0,0],'Linestyle','none','FaceAlpha',0.3); hold on;
    plot(bd2+tt2(zz),att{kk,1}(zz),'r-','Linewidth',2); hold on;
    bb = [(utt{kk,1}(zz)+2*utt{kk,2}(zz)) fliplr(utt{kk,1}(zz)-2*utt{kk,2}(zz))];
    fill(aa,bb,[0,0,1],'Linestyle','none','FaceAlpha',0.3); hold on;
    plot(bd2+tt2(zz),utt{kk,1}(zz),'b-','Linewidth',2); hold on;
    axis tight;
    V = axis;
    axis([-0.05 V(2) 0 vmax]);
    axis off;
    plot([0.005,0.105],[0,0],'k-','Linewidth',2);
    text(0.025,(-0.075*vmax),'100 ms','Fontsize',18);
    plot([0,0],[0.01*vmax,0.75*vmax],'k-');
    text(-0.04,vmax*0.85,'Stim. Onset','Fontsize',16);
    plot([bd,bd],[0.01*vmax,0.75*vmax],'k-');
    text(bd-0.03,vmax*0.85,'Sac. Onset','Fontsize',16);
    plot([-0.04,-0.04],[0,0.2],'k-','Linewidth',2);
    text(-0.05,0.07,'0.2','Rotation',90,'Fontsize',18);
    text(bd+0.005,0.1,sprintf('N=%d',length(zz1)),'Fontsize',16);
    text(-0.08,0.1,'Normed Firing Rate','Rotation',90,'Fontsize',18);
    % text(0.10,-vmax*0.15,'Time','Fontsize',16);
    if (zkk == 1)
       text(0.1,vmax*1.05,'MTC','Fontsize',22,'Color',Ecolo);
    else
       text(0.1,vmax*1.05,'MT','Fontsize',22,'Color',Mcolo); 
    end
  end
  %****************************
  %**** now plot histograms of modulation to the right
  for zkk = 1:3
      % vx = -0.8:0.08:0.8;  %AMY, note here
      vx = -0.5:0.05:0.5;  %AMY, note here
      subplot('Position',[0.65 (0.7-0.3*(zkk-1)) 0.25 0.22]);
      zz1 = find( RATE_Anim == 1);
      zz2 = find( RATE_Anim == 2);
      zzt = find( RATE_Anim > 0);
%      zk = 1;  % Sac Onset resp
      zk = 2;
      if (zkk == 1)
         topname = 'Firing Rate';
         vmax = 30;
         ai1 = (FR{zk}(zz1,3) - FR{zk}(zz1,5) ) ./ (FR{zk}(zz1,3) + FR{zk}(zz1,5));
         ai2 = (FR{zk}(zz2,3) - FR{zk}(zz2,5) ) ./ (FR{zk}(zz2,3) + FR{zk}(zz2,5));                
         ait = (FR{zk}(zzt,3) - FR{zk}(zzt,5) ) ./ (FR{zk}(zzt,3) + FR{zk}(zzt,5));                
      end
      if (zkk == 2)
         topname = 'Fano Factor';
         vmax = 32;
         ai1 = (FF{zk}(zz1,3) - FF{zk}(zz1,5) ) ./ (FF{zk}(zz1,3) + FF{zk}(zz1,5));
         ai2 = (FF{zk}(zz2,3) - FF{zk}(zz2,5) ) ./ (FF{zk}(zz2,3) + FF{zk}(zz2,5)); 
         ait = (FF{zk}(zzt,3) - FF{zk}(zzt,5) ) ./ (FF{zk}(zzt,3) + FF{zk}(zzt,5)); 
      end
      if (zkk == 3)
         topname = 'Mutual Info.';
         vmax = 22;
         ai1 = (MIC{zk}(zz1,3) - MIC{zk}(zz1,5) ) ./ (MIC{zk}(zz1,3) + MIC{zk}(zz1,5));
         ai2 = (MIC{zk}(zz2,3) - MIC{zk}(zz2,5) ) ./ (MIC{zk}(zz2,3) + MIC{zk}(zz2,5));                             
         ait = (MIC{zk}(zzt,3) - MIC{zk}(zzt,5) ) ./ (MIC{zk}(zzt,3) + MIC{zk}(zzt,5));                            
      end
      %*** set white backdrop
      aa = [min(vx) max(vx) max(vx) min(vx) min(vx)];
      bb = [vmax vmax 0 0 vmax];
      fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
      %********
      nyx = hist(ai1,vx);
      nyx = (nyx/sum(nyx))*100;
      byx = hist(ai2,vx);
      byx = (byx/sum(byx))*100;
      bar(vx,nyx,1,'Linestyle','none','FaceColor',Ecolo,'FaceAlpha',0.5); hold on;
      bar(vx,byx,1,'Linestyle','none','FaceColor',Mcolo,'FaceAlpha',0.5);
      axis([min(vx) max(vx) 0 vmax]);
      axis off;
      plot([0,0],[0,vmax],'k-');
      plot([-0.5,-0.5],[0,(0.05*vmax)],'k-','Linewidth',2);
      plot([+0.5,+0.5],[0,(0.05*vmax)],'k-','Linewidth',2);
      plot([min(vx),max(vx)],[0,0],'k-','Linewidth',2);
      if (zkk == 3)
          text(-0.6,-(0.075*vmax),'-0.5','Fontsize',14);
          text(-0.1,-(0.075*vmax),'0.0','Fontsize',14);
          text(0.4,-(0.075*vmax),' 0.5','Fontsize',14);
          %****
          text(-0.3,-(0.20*vmax),'AI Modulation','Fontsize',16);
      end
      plot([min(vx),min(vx)],[1,11],'k-','Linewidth',2);
      text(min(vx)-0.075,4,'10','Rotation',90,'Fontsize',14);
      text(min(vx)-0.2,(0.3*vmax),'Probability','Rotation',90,'Fontsize',16);
      %***
      med1 = nanmedian(ai1); 
      med2 = nanmedian(ai2);
      meda = nanmedian(ait);
      mod1 = 100*(((1+med1)/(1-med1))-1);
      mod2 = 100*(((1+med2)/(1-med2))-1);
      moda = 100*(((1+meda)/(1-meda))-1);
      p1 = signrank(ai1);
      p2 = signrank(ai2);
      p = ranksum(ai1,ai2);
      pa = signrank(ait);
      hg = title(topname);
      set(hg,'Fontsize',16);
      text(max(vx)-0.3,vmax*0.925,'MTC','Fontsize',16,'Color',Ecolo);
      text(max(vx)-0.3,vmax*0.825,'MT','Fontsize',16,'Color',Mcolo);
      % stats in the print out
      % text(max(vx)-0.3,(0.8*vmax),sprintf('E:%5.2f;M:%5.2f',mod1,mod2),'Fontsize',8);
      % text(max(vx)-0.3,(0.7*vmax),sprintf('T:%5.2f(p=%4.3f)',moda,pa),'Fontsize',8);
      disp(topname);
      disp(sprintf('MTC:%5.2f(p=%5.3f);MT:%5.2f(p=%5.3f)',mod1,p1,mod2,p2));
      disp(sprintf('MTC: AI (%5.3f) MT: (%5.3f)',med1,med2));
      disp(sprintf('Total :%5.2f(p=%4.3f)',moda,pa));
      disp(sprintf('Diff between animals p = %10.8f',p));

  end
end
%**********

%*********  make nice figure for paper/thesis
if (1)  % TUNING CURVE FITS
  Mcolo = [0,0.6,0.6];
  Ecolo = [0.3,0.7,0.3];
  hf = figure;
  set(hf,'Position',[100 100 1200 800]);
  %***** first, plot the PSTH per animal on the left
  normo = 1;
  for zkk = 1:2
     %*****
     if (zkk == 1)
      zz1 = find( (RATE_Anim == 1) & (R2List == 1) );  % MTC 1
      subplot('Position',[0.50 0.60 0.20 0.325]);
     else
      zz1 = find( (RATE_Anim == 2) & (R2List == 1) );  % MT 1
      subplot('Position',[0.75 0.60 0.20 0.325]);
     end
     Ibo = VONFIT{1}(zz1,:);  % saccade locked
     %****** get mean curves
     acurves = [];
     ucurves = [];
     maxos = [];
     x = (0:360)*(pi/180);
     xx = 0:360;
     for ck = 1:size(Ibo,1) 
            amu = Ibo(ck,3:6);
            bmu = Ibo(ck,7:10);
            if (1)
              amu(4) = pi;  % fix to standard center
              bmu(4) = pi;  % or could take the average as center...
            else  
              dpref = amu(4)-bmu(4);
              if (dpref > pi)
                dpref = (2*pi)-dpref;
              end
              if (dpref < -pi)
                dpref = dpref + (2*pi);
              end
              amu(4) = pi + (dpref/2);
              bmu(4) = pi - (dpref/2);
            end
            acur = VonMises.vonmises(x,amu); 
            bcur = VonMises.vonmises(x,bmu);
            if (normo) 
              maxo = 0.5*( max(acur) + max(bcur) );
              acur = acur / maxo;  % norm to peak
              bcur = bcur / maxo;  % norm to peak
            end
            acurves = [acurves ; acur];
            ucurves = [ucurves ; bcur];
            maxos = [maxos ; maxo];
     end
     %******
     if (normo)
          mmaxo = 1;
     else
          mmaxo = nanmean(maxos);
     end
     uut = nanmean(ucurves) * mmaxo; 
     suut = (2 * nanstd(ucurves) * mmaxo)/sqrt(length(zz1));
     uut1 = uut + suut;
     uut2 = uut - suut;
     aat = nanmean(acurves) * mmaxo;
     saat = (2 * nanstd(acurves) * mmaxo)/sqrt(length(zz1));
     aat1 = aat + saat;
     aat2 = aat - saat;
     %*****
     vmax = 1.2;
     aa = [0 360 360 0 0];
     bb = [vmax vmax 0 0 vmax];
     fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
     aa = [xx fliplr(xx)];
     bb = [aat1 fliplr(aat2)];
     fill(aa,bb,[1,0,0],'Linestyle','none','FaceAlpha',0.3); hold on;
     plot(xx,aat,'r-','Linewidth',2); hold on;
     bb = [uut1 fliplr(uut2)]; 
     fill(aa,bb,[0,0,1],'Linestyle','none','FaceAlpha',0.3); hold on;
     plot(xx,uut,'b-','Linewidth',2); hold on;
     %******
     axis tight;
     V = axis;
     axis([-15 360 0 vmax]);
     axis off;
     plot([-10,-10],[0,0.2],'k-','Linewidth',2);
     text(-35,0.04,'0.2','Rotation',90,'Fontsize',18);
     text(260,vmax*0.85,sprintf('N=%d',length(zz1)),'Fontsize',16);
     if (zkk == 1)
        text(-70,0.2,'Normed Firing Rate','Rotation',90,'Fontsize',18);
     end
     if (zkk == 1)
         text(215,vmax*0.95,'MTC','Fontsize',18,'Color',Ecolo);
     else
         text(215,vmax*0.95,'MT','Fontsize',18,'Color',Mcolo); 
     end
     plot([90,90],[0,0.05],'k-','Linewidth',2);
     plot([180,180],[0,0.05],'k-','Linewidth',2);
     plot([270,270],[0,0.05],'k-','Linewidth',2);
     plot([0,360],[0,0],'k-','Linewidth',2);
     %******
     text(70,-0.075,'-90','Fontsize',16);
     text(170,-0.075,'0','Fontsize',16);
     text(250,-0.075,'90','Fontsize',16);
     %***
     text(60,-0.20,'Direction (degs)','Fontsize',18);
     %******
  end
  %****************************
  %**** now plot histograms of modulation to the right
  for zkk = 1:3
      vx = -0.8:0.08:0.8; 
      subplot('Position',[0.10 (0.7-0.3*(zkk-1)) 0.30 0.22]);
      zz1 = find( (RATE_Anim == 1) & (R2List == 1) );
      zz2 = find( (RATE_Anim == 2) & (R2List == 1) );
      zzt = find( (RATE_Anim > 0) & (R2List == 1) );
      zk = 1;  % Sac Onset resp
      if (zkk == 1)
         topname = 'Baseline';
         vmax = 30;
         ai1 = (VONBASE{zk}(zz1,3) - VONBASE{zk}(zz1,5) ) ./ (VONBASE{zk}(zz1,3) + VONBASE{zk}(zz1,5));
         ai2 = (VONBASE{zk}(zz2,3) - VONBASE{zk}(zz2,5) ) ./ (VONBASE{zk}(zz2,3) + VONBASE{zk}(zz2,5));                
         ait = (VONBASE{zk}(zzt,3) - VONBASE{zk}(zzt,5) ) ./ (VONBASE{zk}(zzt,3) + VONBASE{zk}(zzt,5));                
      end
      if (zkk == 2)
         topname = 'Gain';
         vmax = 35;
         ai1 = (VONGAIN{zk}(zz1,3) - VONGAIN{zk}(zz1,5) ) ./ (VONGAIN{zk}(zz1,3) + VONGAIN{zk}(zz1,5));
         ai2 = (VONGAIN{zk}(zz2,3) - VONGAIN{zk}(zz2,5) ) ./ (VONGAIN{zk}(zz2,3) + VONGAIN{zk}(zz2,5)); 
         ait = (VONGAIN{zk}(zzt,3) - VONGAIN{zk}(zzt,5) ) ./ (VONGAIN{zk}(zzt,3) + VONGAIN{zk}(zzt,5)); 
      end
      if (zkk == 3)
         topname = 'Width';
         vmax = 40;
         ai1 = (VONKAP{zk}(zz1,3) - VONKAP{zk}(zz1,5) ) ./ (VONKAP{zk}(zz1,3) + VONKAP{zk}(zz1,5));
         ai2 = (VONKAP{zk}(zz2,3) - VONKAP{zk}(zz2,5) ) ./ (VONKAP{zk}(zz2,3) + VONKAP{zk}(zz2,5));                             
         ait = (VONKAP{zk}(zzt,3) - VONKAP{zk}(zzt,5) ) ./ (VONKAP{zk}(zzt,3) + VONKAP{zk}(zzt,5));                            
      end
      %*** set white backdrop
      aa = [min(vx) max(vx) max(vx) min(vx) min(vx)];
      bb = [vmax vmax 0 0 vmax];
      fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
      %********
      nyx = hist(ai1,vx);
      nyx = (nyx/sum(nyx))*100;
      byx = hist(ai2,vx);
      byx = (byx/sum(byx))*100;
      bar(vx,nyx,1,'Linestyle','none','FaceColor',Ecolo,'FaceAlpha',0.5); hold on;
      bar(vx,byx,1,'Linestyle','none','FaceColor',Mcolo,'FaceAlpha',0.5);
      axis([min(vx) max(vx) 0 vmax]);
      axis off;
      plot([0,0],[0,vmax],'k-');
      plot([-0.5,-0.5],[0,(0.05*vmax)],'k-','Linewidth',2);
      plot([+0.5,+0.5],[0,(0.05*vmax)],'k-','Linewidth',2);
      plot([min(vx),max(vx)],[0,0],'k-','Linewidth',2);
      if (zkk == 3)
          text(-0.6,-(0.075*vmax),'-0.5','Fontsize',14);
          text(-0.1,-(0.075*vmax),'0.0','Fontsize',14);
          text(0.4,-(0.075*vmax),' 0.5','Fontsize',14);
          %****
          text(-0.3,-(0.20*vmax),'AI Modulation','Fontsize',18);
      end
      plot([min(vx),min(vx)],[1,11],'k-','Linewidth',2);
      text(min(vx)-0.075,4,'0.1','Rotation',90,'Fontsize',14);
      text(min(vx)-0.2,(0.3*vmax),'Probability','Rotation',90,'Fontsize',18);
      %***
      med1 = nanmedian(ai1); 
      med2 = nanmedian(ai2);
      meda = nanmedian(ait);
      mod1 = 100*(((1+med1)/(1-med1))-1);
      mod2 = 100*(((1+med2)/(1-med2))-1);
      moda = 100*(((1+meda)/(1-meda))-1);
      p = ranksum(ai1,ai2);
      pa = signrank(ait);
      hg = title(topname);
      set(hg,'Fontsize',16);
      text(max(vx)-0.55,vmax*0.925,'MTC','Fontsize',16,'Color',Ecolo);
      text(max(vx)-0.55,vmax*0.825,'MT','Fontsize',16,'Color',Mcolo);
  end
  if (1)
    %*******
    zzt = find( (RATE_Anim == 1) & (R2List == 1) );
    aibase = (VONBASE{zk}(zzt,3) - VONBASE{zk}(zzt,5) ) ./ (VONBASE{zk}(zzt,3) + VONBASE{zk}(zzt,5));                
    aigain = (VONGAIN{zk}(zzt,3) - VONGAIN{zk}(zzt,5) ) ./ (VONGAIN{zk}(zzt,3) + VONGAIN{zk}(zzt,5)); 
    aiwid = (VONKAP{zk}(zzt,3) - VONKAP{zk}(zzt,5) ) ./ (VONKAP{zk}(zzt,3) + VONKAP{zk}(zzt,5));
    aiff = (FF{zk}(zzt,3) - FF{zk}(zzt,5) ) ./ (FF{zk}(zzt,3) + FF{zk}(zzt,5)); 
    aimic = (MIC{zk}(zzt,3) - MIC{zk}(zzt,5) ) ./ (MIC{zk}(zzt,3) + MIC{zk}(zzt,5));                           
    %**** do a linear regression
    XX = [ones(size(aibase)) aibase aigain aiwid aiff];
    [b,bint,r,rint] = regress(aimic,XX);
    aimic2 = XX * b;
    zz = find(~isnan(aimic2));
    br = regress(aimic2(zz),[ones(size(zz)) aimic(zz)]);
    rho = corr(aimic(zz),aimic2(zz));
    %*******
    subplot('Position',[0.50 0.10 0.20 0.325]);
    amax = 0.7;
    aa = [-amax amax amax -amax -amax];
    bb = [-amax -amax amax amax -amax];
    fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
    plot(aimic(zz),aimic2(zz),'k.','Markersize',8); hold on;
    axis off;
    axis([-amax amax -amax amax]);
    % plot([-amax,amax],[-amax,amax],'k-');
    plot([-amax,amax],[(br(1)-br(2)*amax),(br(1)+br(2)*amax)],'k-');
    text(-amax*0.6,amax*0.8,sprintf('R = %5.3f',rho),'Fontsize',14);
    plot([-amax,amax],[-amax,-amax],'k-','Linewidth',2);
    plot([-amax,-amax],[-amax,amax],'k-','Linewidth',2);
    plot([-0.5,-0.5],[-amax,-amax+0.03],'k-','Linewidth',2);
    plot([-0.0,-0.0],[-amax,-amax+0.03],'k-','Linewidth',2);
    plot([+0.5,+0.5],[-amax,-amax+0.03],'k-','Linewidth',2);
    plot([-amax,-amax+0.03],[-0.5,-0.5],'k-','Linewidth',2);
    plot([-amax,-amax+0.03],[-0.0,-0.0],'k-','Linewidth',2);
    plot([-amax,-amax+0.03],[+0.5,+0.5],'k-','Linewidth',2);
    text(-0.575,-amax-0.1,'-0.5','Fontsize',14);
    text(-0.05,-amax-0.1,' 0.0','Fontsize',14);
    text(+0.425,-amax-0.1,'+0.5','Fontsize',14);
    text(-amax-0.1,-0.575,'-0.5','Fontsize',14,'Rotation',90);
    text(-amax-0.1,-0.05,' 0.0','Fontsize',14,'Rotation',90);
    text(-amax-0.1,0.425,'+0.5','Fontsize',14,'Rotation',90);
    text(-0.4,-amax-0.25,'AI Mutual Info.','Fontsize',18);
    text(-amax-0.3,-0.4,'Linear Prediction','Fontsize',18,'Rotation',90);
    %********
    %********
    subplot('Position',[0.80 0.10 0.175 0.325]);
    aa = [1 6 6 1 1];
    bb = [-0.6 -0.6 0.9 0.9 -0.6];
    fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
    bar(2:5,b(2:5),1,'Linestyle','none','FaceColor',[0.5,0.5,0.5]); hold on;
    for zk = 2:length(b)
       plot([zk,zk],[bint(zk,1),bint(zk,2)],'k-','Linewidth',2); 
    end
    axis([1 6 -0.6 0.9]);
    axis off;
    plot([1.2,5.7],[0,0],'k-','Linewidth',1);
    plot([1.1,1.1],[0,0.5],'k-','Linewidth',2);
    text(0.8,0.2,'0.5','Rotation',90,'Fontsize',14);
    text(2,-0.6,'Base','Rotation',90,'Fontsize',16);
    text(3,-0.6,'Gain','Rotation',90,'Fontsize',16);
    text(4,-0.6,'Width','Rotation',90,'Fontsize',16);
    text(5,-0.6,'FF','Rotation',90,'Fontsize',16);
    text(0.1,-0.3,'Linear Coefficient','Rotation',90,'Fontsize',18);
    %*******  
  end
end
%**********

%*********  make nice figure for paper/thesis
if (1)  % NON-PARAMETRIC STATS FOR TUNING CURVE FITS
  Mcolo = [0,0.6,0.6];
  Ecolo = [0.3,0.7,0.3];
  hf = figure;
  set(hf,'Position',[100 50 900 900]);
  %****************************
  if (1)  % First show that von mises modulations are correlated (possible biases?)
    %*******
    zk = 1;
    zzt = find( (RATE_Anim > 0) & (R2List == 1));  % only good R2 fits
    aibase = (VONBASE{zk}(zzt,3) - VONBASE{zk}(zzt,5) ) ./ (VONBASE{zk}(zzt,3) + VONBASE{zk}(zzt,5));                
    aigain = (VONGAIN{zk}(zzt,3) - VONGAIN{zk}(zzt,5) ) ./ (VONGAIN{zk}(zzt,3) + VONGAIN{zk}(zzt,5)); 
    aiwid = (VONKAP{zk}(zzt,3) - VONKAP{zk}(zzt,5) ) ./ (VONKAP{zk}(zzt,3) + VONKAP{zk}(zzt,5));
    
    %aiwid = 0.5*(VONKAP{zk}(zzt,3) + VONKAP{zk}(zzt,5) )/180;
    %********
    nzz = find( (NPBASE(:,1) > 0) ); 
    naibase = (NPBASE(nzz,2) - NPBASE(nzz,4) ) ./ (NPBASE(nzz,2) + NPBASE(nzz,4));                
    naigain = (NPAMP(nzz,2) - NPAMP(nzz,4) ) ./ (NPAMP(nzz,2) + NPAMP(nzz,4)); 
    naiwid = (NPWID(nzz,2) - NPWID(nzz,4) ) ./ (NPWID(nzz,2) + NPWID(nzz,4));
    nailoc = LOCBIAS(nzz,1);
    
    if (0) % factor out ailoc
      b = regress(naibase,nailoc);
      yhat = b(1) * nailoc;
      naibase = naibase - yhat;  % regress out factor
    end
    %*******
    rnzz = find( (RNPBASE(:,1)>0) ); % & (RNPBASE(:,2) > 0) & (RNPBASE(:,4) > 0) & ...
                  % (RNPAMP(:,2) > RNPAMP(:,4)) );
    rnaibase = (RNPBASE(rnzz,2) - RNPBASE(rnzz,4) ) ./ (RNPBASE(rnzz,2) + RNPBASE(rnzz,4));                
    rnaigain = (RNPAMP(rnzz,2) - RNPAMP(rnzz,4) ) ./ (RNPAMP(rnzz,2) + RNPAMP(rnzz,4)); 
    rnaiwid = (RNPWID(rnzz,2) - RNPWID(rnzz,4) ) ./ (RNPWID(rnzz,2) + RNPWID(rnzz,4));
    %******
    aiff = (FF{zk}(nzz,3) - FF{zk}(nzz,5) ) ./ (FF{zk}(nzz,3) + FF{zk}(nzz,5)); 
    aimic = (MIC{zk}(nzz,3) - MIC{zk}(nzz,5) ) ./ (MIC{zk}(nzz,3) + MIC{zk}(nzz,5));                           
    
    %****** NON-PARAMETRIC STATS
    zt1 = find(RATE_Anim(nzz) == 1);
    zt2 = find(RATE_Anim(nzz) == 2);
    disp('  ');
    disp('  ');
    disp('Non parametric fit AI indices');
    disp(sprintf('MTC ( N = %d)  MT (N = %d)',length(naibase(zt1)),length(naibase(zt2))));
    compute_ai_print(naibase(zt1),naibase(zt2),'Baseline rates ');
    compute_ai_print(naigain(zt1),naigain(zt2),'Amplitudes ');
    compute_ai_print(naiwid(zt1),naiwid(zt2),'Widths ');
    %**** USE THE MORE RESTRICTED VON MISES SET (zzt instead of nzz)
    vt1 = find(RATE_Anim(zzt) == 1);
    vt2 = find(RATE_Anim(zzt) == 2);
    disp('  ');
    disp('  ');
    disp('Von mises fits AI indices');
    disp(sprintf('MTC ( N = %d)  MT (N = %d)',length(aibase(vt1)),length(aibase(vt2))));
    compute_ai_print(aibase(vt1),aibase(vt2),'Baseline rates ');
    compute_ai_print(aigain(vt1),aigain(vt2),'Amplitudes ');
    compute_ai_print(aiwid(vt1),aiwid(vt2),'Widths ');
    
    input('check stats');
    
    disp('  ');
    if (0)
      rzt1 = find(RATE_Anim(rnzz) == 1);
      rzt2 = find(RATE_Anim(rnzz) == 2);
      disp(' Null NP random shuffles');
      disp('   ');
      compute_ai_print(rnaibase(rzt1),rnaibase(rzt2),'R Baseline rates ');
      compute_ai_print(rnaigain(rzt1),rnaigain(rzt2),'R Amplitudes ');
      compute_ai_print(rnaiwid(rzt1),rnaiwid(rzt2),'R Widths ');
      disp('   ');
    end
    %**** do a linear regression
    amax = 1.0;
    for ik = 1:9  % three scatter plots to look for von mises biases
        %*********
        if (ik == 1)
           % subplot('Position',[0.08 0.725 0.225 0.225]);
           xname = 'AI Baseline';
           yname = 'AI Gain';
           zz = find( ~isnan(aibase) & ~isnan(aigain) );
           ai1 = aibase(zz);
           ai2 = aigain(zz);
        end
        if (ik == 2)
           % subplot('Position',[0.40 0.725 0.225 0.225]);
           xname = 'AI Gain';
           yname = 'AI Width';
           zz = find( ~isnan(aigain) & ~isnan(aiwid) );
           ai1 = aigain(zz);
           ai2 = aiwid(zz);
        end
        if (ik == 3)
           % subplot('Position',[0.72 0.725 0.225 0.225]);
           xname = 'AI Width';
           yname = 'AI Baseline';
           zz = find( ~isnan(aiwid) & ~isnan(aibase) );
           ai1 = aiwid(zz);
           ai2 = aibase(zz);
        end
        if (ik == 7)
           % subplot('Position',[0.08 0.725 0.225 0.225]);
           xname = 'AI Baseline';
           yname = 'AI Location';
           zz = find( ~isnan(naibase) & ~isnan(nailoc) ); % & (NPBASE(nzz,1)==2) );
           ai1 = naibase(zz);
           ai2 = nailoc(zz);
        end
        if (ik == 8)
           % subplot('Position',[0.40 0.725 0.225 0.225]);
           xname = 'AI Gain';
           yname = 'AI Location';
           zz = find( ~isnan(naigain) & ~isnan(nailoc) ); % & (NPBASE(nzz,1)==2) );
           ai1 = naigain(zz);
           ai2 = nailoc(zz);
        end
        if (ik == 9)
           % subplot('Position',[0.72 0.725 0.225 0.225]);
           xname = 'AI Width';
           yname = 'AI Location';
           zz = find( ~isnan(naiwid) & ~isnan(nailoc) ); % & (NPBASE(nzz,1)==2) );
           ai1 = naiwid(zz);
           ai2 = nailoc(zz);
        end
        %********* NP parameters *******
        if (ik == 4)
           subplot('Position',[0.08 0.33 0.225 0.225]);
           xname = 'AI Baseline';
           yname = 'AI Gain';
           xoff = 0.15;
           zz = find( ~isnan(naibase) & ~isnan(naigain) ); % & (NPBASE(nzz,1)==2) );
           ai1 = naibase(zz);
           ai2 = naigain(zz);
        end
        if (ik == 5)
           subplot('Position',[0.40 0.33 0.225 0.225]);
           xname = 'AI Gain';
           yname = 'AI Width';
           xoff = 0.2;
           zz = find( ~isnan(naigain) & ~isnan(naiwid) );
           ai1 = naigain(zz);
           ai2 = naiwid(zz);
        end
        if (ik == 6)
           subplot('Position',[0.72 0.33 0.225 0.225]);
           xname = 'AI Width';
           yname = 'AI Baseline';
           xoff = 0.3;
           zz = find( ~isnan(naiwid) & ~isnan(naibase) );
           ai1 = naiwid(zz);
           ai2 = naibase(zz);
        end
        %*********
        [rho,pho] = corr(ai1,ai2,'type','Spearman');  % non-parametric
        if (ik > 3)
            if (ik < 7)
               disp(sprintf('NP %s x %s r=%5.3f (p=%5.3f)',xname,yname,rho,pho));
            else
               disp(sprintf('CONTROL: **** NP %s x AI LOCATION r=%5.3f (p=%5.3f)',xname,rho,pho));    
            end
        else
            disp(sprintf('Von Mises %s x %s r=%5.3f (p=%5.3f)',xname,yname,rho,pho));    
        end
        %****
        if (ik > 3) && (ik < 7)
          aa = [-amax amax amax -amax -amax];
          bb = [-amax -amax amax amax -amax];
          fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
          zta = find( RATE_Anim(zz) == 2);
          plot(ai1(zta),ai2(zta),'k.','Markersize',8,'Color',Mcolo); hold on;
          zta = find( RATE_Anim(zz) == 1);
          plot(ai1(zta),ai2(zta),'k.','Markersize',8,'Color',Ecolo); hold on;
          % plot(ai1,ai2,'k.','Markersize',8); hold on;
          if (1)
           plot([-amax,amax],[0,0],'k-','Linewidth',1,'Color',[0.7,0.7,0.7]);
           plot([0,0],[-amax,amax],'k-','Linewidth',1,'Color',[0.7,0.7,0.7]);           
          end
          axis off;
          axis([-amax amax -amax amax]);
          text(amax*0.25,amax*1.22,sprintf('R = %5.3f',rho),'Fontsize',12);
          text(amax*0.25,amax*1.08,sprintf('  p = %5.3f',pho),'Fontsize',12);
          plot([-amax,amax],[-amax,-amax],'k-','Linewidth',1);
          plot([-amax,-amax],[-amax,amax],'k-','Linewidth',1);
          plot([-0.5,-0.5],[-amax,-amax+0.06],'k-','Linewidth',2);
          plot([-0.0,-0.0],[-amax,-amax+0.06],'k-','Linewidth',2);
          plot([+0.5,+0.5],[-amax,-amax+0.06],'k-','Linewidth',2);
          plot([-amax,-amax+0.06],[-0.5,-0.5],'k-','Linewidth',2);
          plot([-amax,-amax+0.06],[-0.0,-0.0],'k-','Linewidth',2);
          plot([-amax,-amax+0.06],[+0.5,+0.5],'k-','Linewidth',2);
          % text(-0.7,-amax-0.1,'-0.5','Fontsize',14);
          % text(-0.2,-amax-0.1,' 0.0','Fontsize',14);
          % text(+0.3,-amax-0.1,'+0.5','Fontsize',14);
          text(-amax-0.1,-0.75,'-0.5','Fontsize',14,'Rotation',90);
          text(-amax-0.1,-0.25,' 0.0','Fontsize',14,'Rotation',90);
          text(-amax-0.1,0.25,'+0.5','Fontsize',14,'Rotation',90);
          % text(-0.4,amax*1.2,xname,'Fontsize',16);
          text(-amax-0.4,-0.2-xoff,yname,'Fontsize',16,'Rotation',90);
          if (ik == 2)
             % text(-amax,amax*1.2,'Von Mises Parameters','Fontsize',18);
          end
          if (ik == 5)
             % text(-1.5*amax,amax*1.4,'Non-Parametric Tuning Statistics','Fontsize',18);         
          end
        end
    end
  end  

  
  
  
  %**** now plot histograms of modulation to the left
  for zkk = 1:3
      vx = -amax:(amax/10):amax; 
      subplot('Position',[(0.08+(0.32*(zkk-1))) 0.17 0.225 0.15]);
      % subplot('Position',[0.08 (0.44-((zkk-1)*0.19)) 0.225 0.15]);
      zz1 = find( (RATE_Anim == 1) ); % & ...
                       % (NPBASE(:,2)>0) & (NPBASE(:,4)>0) & ...
                       % (NPAMP(:,2) > NPAMP(:,4)) );  %& (R2List == 1) );
      zz2 = find( (RATE_Anim == 2) ); % & ...
                       % (NPBASE(:,2)>0) & (NPBASE(:,4)>0) & ...
                       % (NPAMP(:,2) > NPAMP(:,4)) ); %& (R2List == 1) );
      zzt = find( (RATE_Anim > 0) ); % & ...
                       % (NPBASE(:,2)>0) & (NPBASE(:,4)>0) & ...
                       % (NPAMP(:,2) > NPAMP(:,4)) ); %& (R2List == 1) );
      if (zkk == 1)
         topname = 'AI Baseline';
         xoff = 0.0;
         vmax = 25;
         ai1 = (NPBASE(zz1,2) - NPBASE(zz1,4) ) ./ (NPBASE(zz1,2) + NPBASE(zz1,4));
         ai2 = (NPBASE(zz2,2) - NPBASE(zz2,4) ) ./ (NPBASE(zz2,2) + NPBASE(zz2,4));                
         ait = (NPBASE(zzt,2) - NPBASE(zzt,4) ) ./ (NPBASE(zzt,2) + NPBASE(zzt,4));                
      end
      if (zkk == 2)
         topname = 'AI Gain';
         xoff = -0.2;
         vmax = 30;
         ai1 = (NPAMP(zz1,2) - NPAMP(zz1,4) ) ./ (NPAMP(zz1,2) + NPAMP(zz1,4));
         ai2 = (NPAMP(zz2,2) - NPAMP(zz2,4) ) ./ (NPAMP(zz2,2) + NPAMP(zz2,4)); 
         ait = (NPAMP(zzt,2) - NPAMP(zzt,4) ) ./ (NPAMP(zzt,2) + NPAMP(zzt,4)); 
      end
      if (zkk == 3)
         topname = 'AI Width';
         xoff = -0.1;
         vmax = 35;
         ai1 = (NPWID(zz1,2) - NPWID(zz1,4) ) ./ (NPWID(zz1,2) + NPWID(zz1,4));
         ai2 = (NPWID(zz2,2) - NPWID(zz2,4) ) ./ (NPWID(zz2,2) + NPWID(zz2,4));                             
         ait = (NPWID(zzt,2) - NPWID(zzt,4) ) ./ (NPWID(zzt,2) + NPWID(zzt,4));                            
      end
      %*** set white backdrop
      aa = [min(vx) max(vx) max(vx) min(vx) min(vx)];
      bb = [vmax vmax 0 0 vmax];
      fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
      %********
      nyx = hist(ai1,vx);
      nyx = (nyx/sum(nyx))*100;
      byx = hist(ai2,vx);
      byx = (byx/sum(byx))*100;
      bar(vx,nyx,1,'Linestyle','none','FaceColor',Ecolo,'FaceAlpha',0.5); hold on;
      bar(vx,byx,1,'Linestyle','none','FaceColor',Mcolo,'FaceAlpha',0.5);
      axis([min(vx) max(vx) 0 vmax]);
      axis off;
      plot([0,0],[0,vmax],'k-');
      plot([-0.5,-0.5],[0,(0.05*vmax)],'k-','Linewidth',2);
      plot([+0.5,+0.5],[0,(0.05*vmax)],'k-','Linewidth',2);
      plot([min(vx),max(vx)],[0,0],'k-','Linewidth',2);
      if (1) % (zkk == 3)
          text(-0.70,-(0.075*vmax),'-0.5','Fontsize',14);
          text(-0.10,-(0.075*vmax),'0.0','Fontsize',14);
          text(0.35,-(0.075*vmax),' 0.5','Fontsize',14);
          %****
          text(-0.55-xoff,-(0.25*vmax),topname,'Fontsize',16); % 'AI Modulation','Fontsize',16);
      end
      plot([min(vx),min(vx)],[1,11],'k-','Linewidth',2);
      text(min(vx)-0.125,3,'0.1','Rotation',90,'Fontsize',14);
      text(min(vx)-0.4,(0.15*vmax),'Probability','Rotation',90,'Fontsize',16);
      %***
      med1 = nanmedian(ai1); 
      med2 = nanmedian(ai2);
      meda = nanmedian(ait);
      mod1 = 100*(((1+med1)/(1-med1))-1);
      mod2 = 100*(((1+med2)/(1-med2))-1);
      moda = 100*(((1+meda)/(1-meda))-1);
      p = ranksum(ai1,ai2);
      pa = signrank(ait);
      % hg = title(topname);
      % set(hg,'Fontsize',16);
      text(max(vx)-0.65,vmax*0.90,'MTC','Fontsize',12,'Color',Ecolo);
      text(max(vx)-0.65,vmax*0.80,'MT','Fontsize',12,'Color',Mcolo);
  end
  %******** last would plot a single unit example non-parametric curve
  %******* plot the average over population
  if (1)
     %******* example of NP computation for unit in Fig 3 
     %************
     %Info = ExampleInfo;
     afit1 = Info.SacOn.OriTune{1};
     ovals = Info.SacOn.OriTune{1}.OriVals;
     ovals = (ovals+22.5)/22.5;
     afit2 = Info.SacOn.OriTune{2};
     LO = length(ovals);
     %****** find prefs
     % smooth on the circle ******
     aa = circsmooth(afit1.Mu,NPSMOOTH);
     uu = circsmooth(afit2.Mu,NPSMOOTH);
     sa = circsmooth(afit1.Sem,NPSMOOTH);
     su = circsmooth(afit2.Sem,NPSMOOTH);
     tot = aa + uu;
     theta = 0:(LO-1);
     theta = theta*(2*pi)/LO;
     wvec = sum( tot .* exp(j*theta) );
     wvec = wvec/abs(wvec);
     %********
     if (centerwid) 
          awvec = sum( aa .* exp(j*theta) );
          awvec = awvec/abs(awvec);
          uwvec = sum( uu .* exp(j*theta) );
          uwvec = uwvec/abs(uwvec);
          %*******
          wvec = 0.5*(awvec+uwvec);
          wvec = wvec / abs(wvec);
     end
     %********
     Wvec2 = wvec;
     ango = angle(wvec);
     if (ango < 0)
         ango = ango + 2*pi;
     end
     wvc = (LO*ango/(2*pi)) + 1;
     wshift = 3;
     wvc = wvc + wshift;
     rval = [];
     for ii = 1:length(oval1)
                ang = ((ii-1)*2*pi/length(oval1));
                dot = (real(Wvec2)*cos(ang) + imag(Wvec2)*sin(ang));
                rval = [rval ; [ii dot 1 cos(ang) sin(ang)]];
     end
     yrvals = sortrows(rval,2);
     rord = yrvals(:,1);
     NPSMOOTH = 1;
     ratt = aa(rord);  %smoothnp(afit1.Mu(rord),NPSMOOTH);
     sratt = sa(rord); %smoothnp(afit1.Sem(rord),NPSMOOTH);
     rutt = uu(rord);  %smoothnp(afit2.Mu(rord),NPSMOOTH);
     srutt = su(rord); %smoothnp(afit2.Sem(rord),NPSMOOTH);
     [ampw,umpw,abasew,ubasew,aw,uw] = compute_NP_params(ratt,rutt);
     %***********       
     vmax = 1.1*max([max(afit1.Mu),max(afit2.Mu)]);
     aa = [0 (LO+1) (LO+1) 0 0];
     bb = [vmax vmax 0 0 vmax];
     %******* plot original tuning
     subplot('Position',[0.08 0.810 0.225 0.125]);
     fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
     axis([0 LO+1 0 vmax]);
     axis off;
     xx = 1:LO;
     if (wshift > 0)
        xx = [(LO-wshift+1):LO, 1:(LO-wshift)]; % recenter
     end
     %******
     aa = circsmooth(afit1.Mu,NPSMOOTH);
     uu = circsmooth(afit2.Mu,NPSMOOTH);
     sa = circsmooth(afit1.Sem,NPSMOOTH);
     su = circsmooth(afit2.Sem,NPSMOOTH);
     %******
     zaa = [1:LO fliplr(1:LO)];
     zbb = [(aa(xx)+2*sa(xx)) fliplr(aa(xx)-2*sa(xx))];
     fill(zaa,zbb,[1,0,0],'FaceAlpha',0.1,'Linestyle','none'); hold on;
     zaa = [1:LO fliplr(1:LO)];
     zbb = [(uu(xx)+2*su(xx)) fliplr(uu(xx)-2*su(xx))];
     fill(zaa,zbb,[0,0,1],'FaceAlpha',0.1,'Linestyle','none'); hold on;
     plot(ovals,aa(xx),'r.-','Markersize',10); hold on;
     plot(ovals,uu(xx),'b.-','Markersize',10); hold on;
     plot([wvc,wvc],[0,vmax],'k--','Linewidth',1);
     text(4,-0.1*vmax,'Motion direction','Fontsize',14);
     plot([0.0,0.0],[0,10],'k-','Linewidth',2);
     text(-1.0,3,'10','Fontsize',12,'Rotation',90);
     text(-2.5,0.1*vmax,'Rate (sp/s)','Rotation',90,'Fontsize',14);
     %*****
     text(3,1.40*vmax,'Example Unit','Fontsize',16);
     text(-2,1.20*vmax,'Non-Parametric(NP) Tuning','Fontsize',16);
     plot([1,LO],[0,0],'k-');
     
     %****** plot tuning ranked by preference ****
     subplot('Position',[0.08 0.65 0.225 0.125]); 
     aa = [-1 (LO+2) (LO+2) -1 -1];
     fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
     axis([-1 (LO+2) 0 vmax]);
     axis off;
     if (1)
       plot([15.0,15.0],[0,abasew],'r-','LineWidth',2,'Color',[1,0,0]); hold on;
       plot([15.5,15.5],[0,ubasew],'b-','LineWidth',2,'Color',[0,0,1]);
       plot([17,17],[abasew,abasew+ampw],'r-','LineWidth',2,'Color',[1,0,0]);
       plot([17.5,17.5],[ubasew,ubasew+umpw],'b-','LineWidth',2,'Color',[0,0,1]);
       plot([16-aw,16],[0.975*vmax,0.975*vmax],'r-','LineWidth',2,'Color',[1,0,0]);
       plot([16-uw,16],[0.925*vmax,0.925*vmax],'b-','LineWidth',2,'Color',[0,0,1]);
     end
     plot([-2,(LO+2)],[0,0],'k-');
     %*****
     aa = [1:LO fliplr(1:LO)];
     bb = [(ratt+2*sratt) fliplr(ratt-2*sratt)];
     fill(aa,bb,[1,0,0],'FaceAlpha',0.1,'Linestyle','none'); hold on;
     plot(1:LO,ratt,'r.-','Markersize',10); hold on;
     aa = [1:LO fliplr(1:LO)];
     bb = [(rutt+2*srutt) fliplr(rutt-2*srutt)];
     fill(aa,bb,[0,0,1],'FaceAlpha',0.1,'Linestyle','none'); hold on;
     plot(1:LO,rutt,'b.-','Markersize',10);
     %*******
     wvv = min([wvc-floor(wvc),ceil(wvc)-wvc]);
     plot([LO+wvv,LO+wvv],[0,vmax],'k--','Linewidth',1);
     text(-1,-0.125*vmax,'Direction by Preference','Fontsize',14);
     plot([-1,-1],[0,10],'k-','Linewidth',2);
     text(-2,3,'10','Fontsize',12,'Rotation',90);
     text(-4.5,0.1*vmax,'Rate (sp/s)','Rotation',90,'Fontsize',14);
     
     %****************
  end   
  if (1)
     subplot('Position',[0.40 0.65 0.225 0.30]);
     %**** get average tuning shapes
     aarate = [];
     uurate = [];
     for zk = 1:2
       if (zk == 1)
           zz = find( (RATE_Anim == 1) ); % & ...
                       % (NPBASE(:,2)>0) & (NPBASE(:,4)>0) & ...
                       % (NPAMP(:,2) > NPAMP(:,4)) );
       end
       if (zk == 2)
           zz = find( (RATE_Anim == 2) ); % & ...
                       % (NPBASE(:,2)>0) & (NPBASE(:,4)>0) & ...
                       % (NPAMP(:,2) > NPAMP(:,4)) );
       end
       %*****
       arate = [];
       irate = [];
       for znp = 1:size(zz,1) % NPFIT{1},1)
          ma = [];
          mb = [];
          for uk = 1:2  % do over both conditions so we can keep on same scale
             ma = [ma max(NPFIT{uk}(zz(znp),:))];
          end
          maxo = mean(ma);
          aa = ( NPFIT{1}(zz(znp),:)/maxo);
          uu = ( NPFIT{2}(zz(znp),:)/maxo);
          %**********
          arate = [arate ; aa];
          irate = [irate ; uu];
       end
       aarate = [aarate ; arate];
       uurate = [uurate ; irate];
       %*****
       TT{zk,1} = 1:LO;
       TT{zk,2} = TT{zk,1};
       UU{zk,1} = nanmean(arate);
       SU{zk,1} = nanstd(arate)/sqrt(size(arate,1));
       UU{zk,2} = nanmean(irate);
       SU{zk,2} = nanstd(irate)/sqrt(size(irate,1));
     end
     %******* weight two monkeys equal in average
     if (0)
       aTT{1} = 0.5*(TT{1,1}+TT{2,1});
       aTT{2} = 0.5*(TT{1,2}+TT{2,2});
       aUU{1} = 0.5*(UU{1,1}+UU{2,1});
       aUU{2} = 0.5*(UU{1,2}+UU{2,2});
       aSU{1} = 0.5*(SU{1,1}+SU{2,1});
       aSU{2} = 0.5*(SU{1,2}+SU{2,2});
     else
       aTT{1} = 0.5*(TT{1,1}+TT{2,1});
       aTT{2} = 0.5*(TT{1,2}+TT{2,2});
       aUU{1} = nanmean(aarate);
       aSU{1} = nanstd(aarate)/sqrt(size(aarate,1));
       aUU{2} = nanmean(uurate);
       aSU{2} = nanstd(uurate)/sqrt(size(uurate,1));
     end
     %********
     vmax = 1;
     LO = size(arate,2);
     aa = [0 LO LO 0 0];
     bb = [vmax vmax 0 0 vmax];
     fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
     axis([0 LO 0 vmax]);
     axis off;
     %******
     text(2,1.1*vmax,'Mean NP Tuning','Fontsize',16);         
     %*******
     spikeplot.plot_attention_traces(aTT,aUU,aSU,[1,LO],[],[]);
     axis([0 LO 0 vmax]);
     axis off;
     plot([(LO-7),LO],[0,0],'k-','Linewidth',2);
     text((LO-6),0.1*vmax,'180 degs','Fontsize',14);
     plot([0.5,0.5],[0,0.5],'k-','Linewidth',2);
     text(-0.4,0.2,'0.5','Rotation',90,'Fontsize',14);
     text(-3.0,0.2,'Normalized Rate','Rotation',90,'Fontsize',14);
     text(0,-0.075*vmax,'Direction by Preference','Fontsize',14);
     text(2,0.85*vmax,sprintf('N=%d',length(aarate)),'Fontsize',14);
     %*******
  end
  %******** and plot -log pvalues as histograms
  if (0)
     %******* plot distribution of t-values
     if (1) % use NP stats here
         zz = find(~isnan(NPWID(:,2)) & ~isnan(NPWID(:,4)) );
         zz1 = find( (NPWID(zz,1) == 1) ); % & ...
                     % (NPBASE(zz,2)>0) & (NPBASE(zz,4)>0) ); % MTC
         zz2 = find( (NPWID(zz,1) == 2) ); % & ...
                     % (NPBASE(zz,2)>0) & (NPBASE(zz,4)>0) ); % MT
         comp1 = abs( NPWID(zz(zz1),2) - NPWID(zz(zz1),4) );
         comp2 = abs( NPWID(zz(zz2),2) - NPWID(zz(zz2),4) );
         wido1 = sqrt( ( NPWID(zz(zz1),3).^2 + NPWID(zz(zz1),5).^2));
         wido2 = sqrt( ( NPWID(zz(zz2),3).^2 + NPWID(zz(zz2),5).^2));
     else
         %aiwid = (VONKAP{zk}(zzt,3) - VONKAP{zk}(zzt,5) ) ./ (VONKAP{zk}(zzt,3) + VONKAP{zk}(zzt,5));
 
         zz = find(R2List == 1);
         zk = 1;
         zz1 = find( (VONKAP{zk}(zz,1) == 1) );
         zz2 = find( (VONKAP{zk}(zz,1) == 2) );  
         comp1 = abs( VONKAP{zk}(zz(zz1),3) - VONKAP{zk}(zz(zz1),5) );
         comp2 = abs( VONKAP{zk}(zz(zz2),3) - VONKAP{zk}(zz(zz2),5) );
         wido1 = sqrt( ( VONKAP{zk}(zz(zz1),4).^2 + VONKAP{zk}(zz(zz1),6).^2));
         wido2 = sqrt( ( VONKAP{zk}(zz(zz2),4).^2 + VONKAP{zk}(zz(zz2),6).^2));
    
     end
     tval1 = comp1 ./ wido1;
     tval2 = comp2 ./ wido2;
     tcrit = 1.96; % 2.0;
     ai1 = (NPWID(zz(zz1),2)-NPWID(zz(zz1),4))./(NPWID(zz(zz1),2)+NPWID(zz(zz1),4));
     ai2 = (NPWID(zz(zz2),2)-NPWID(zz(zz2),4))./(NPWID(zz(zz2),2)+NPWID(zz(zz2),4));
     %******** random shuffle (null distribution)
     zz = find( ~isnan(RNPWID(:,2)) & ~isnan(RNPWID(:,4)) );
     rcomp = abs( RNPWID(zz,2) - RNPWID(zz,4) );
     rwido = sqrt( ( RNPWID(zz,3).^2 + RNPWID(zz,5).^2));
     rtval = rcomp ./ rwido;
     
     %***** get stats on how likely to see this many sig 
     P1 = sum(tval1 > tcrit);
     M1 = length(zz1);
     pval1 = binocdf(P1,M1,0.05,'upper');
     disp(sprintf('Marm E, N=%d, M=%d, perc %4.2f p = %5.4f',M1,P1,(100*P1/M1),pval1));
     P2 = sum(tval2 > tcrit);
     M2 = length(zz2);
     pval2 = binocdf(P2,M2,0.05,'upper');
     disp(sprintf('Marm M, N=%d, M=%d, perc %4.2f p = %5.4f',M2,P2,(100*P2/M2),pval2));
     disp(sprintf('Total percentatge %4.2f',(100*(P1+P2)/(M1+M2))));
     
     %*******
     subplot('Position',[0.72 0.65 0.225 0.30]);
     pmax = 40;
     xmax = 8;
     aa = [0 xmax xmax 0 0];
     bb = [0 0 pmax pmax 0];
     fill(aa,bb,[1,1,1],'Linestyle','none'); hold on;
     %******
     vx = 0:0.4:xmax;
     yx1 = hist(tval1,vx);
     yx1 = 100*yx1/sum(yx1);
     yx2 = hist(tval2,vx);
     yx2 = 100*yx2/sum(yx2);
     yxt = hist(rtval,vx);
     yxt = 100*yxt/sum(yxt);
     %*****
     bar(vx,yx1,1,'Linestyle','none','FaceColor',Ecolo,'FaceAlpha',0.5); hold on;
     bar(vx,yx2,1,'Linestyle','none','FaceColor',Mcolo,'FaceAlpha',0.5);
     plot(vx,yxt,'k-','Linewidth',1);
     plot([tcrit,tcrit],[0,pmax],'k--');
     text(tcrit-0.3,-0.05*pmax,'2','Fontsize',14);
     axis([0 xmax 0 pmax]);
     axis off;
     text(0.4*xmax,-0.075*pmax,'Z-value','Fontsize',14);
     plot([0,0],[0,10],'k-','Linewidth',2);
     text(-1,3,'10','Fontsize',14,'Rotation',90);
     text(-2,0.3*pmax,'Probability','Rotation',90,'Fontsize',14);
     text(0.0*xmax,1.1*pmax,'Sig. Width Changes','Fontsize',16);
     sig1 = 100*(sum(tval1 > tcrit)/length(tval1));
     text(0.35*xmax,0.9*pmax,[sprintf('MTC %4.1f',sig1),'%'],...
                          'Fontsize',14,'Color',Ecolo);
     sig2 = 100*(sum(tval2 > tcrit)/length(tval2));
     text(0.35*xmax,0.82*pmax,[sprintf('Marm M %4.1f',sig2),'%'],...
                         'Fontsize',14,'Color',Mcolo);
     % sigt = 100*(sum(rtval > tcrit)/length(rtval));
     % text(0.35*xmax,0.74*pmax,[sprintf('Null Dist. %4.1f',sigt),'%'],...
     %                    'Fontsize',14,'Color',[0.3,0.3,0.3]);
                  
     plot([0,xmax],[0,0],'k-');
     %*******     
  end
end

showall = 1;
if (0) || (showall == 1)   %****** now plot the psth plots ******
  normo = 1;
  Anim = VONFIT{1}(:,1); % MTC, MT, or Laminar? 
  hf = figure(22);
  set(hf,'Position',[100 150 1200 800]);
  TT = cell(3,2);  UU = cell(3,2);  SU = cell(3,2); % 1-3 animal, 1-2 timeint
  %*************
  for nn = 1:3  
    zz1 = find( RATE_Anim == nn);  % MTC 1, MT Tung 2, MT Laminar 3
    for kk = 1:2 
      ha = subplot('Position',[(0.1+(0.3*(nn-1))) (0.1+(0.5*(kk-1))) 0.20 0.35]);
      if nn == 3
          %****** pool the results from the previous two and plot them
          TT{3,kk}{1} = 0.5*(TT{1,kk}{1} + TT{2,kk}{1});
          TT{3,kk}{2} = 0.5*(TT{1,kk}{2} + TT{2,kk}{2});
          UU{3,kk}{1} = 0.5*(UU{1,kk}{1} + UU{2,kk}{1});
          SU{3,kk}{1} = 0.5*(SU{1,kk}{1} + SU{2,kk}{1});
          UU{3,kk}{2} = 0.5*(UU{1,kk}{2} + UU{2,kk}{2});
          SU{3,kk}{2} = 0.5*(SU{1,kk}{2} + SU{2,kk}{2});
          zz1 = [1,2];  % two animals
          %******        
      else
          TT{nn,kk}{1} = nanmean(RATE_tt{kk}(zz1,:));
          TT{nn,kk}{2} = TT{nn,kk}{1};
          if (normo == 1) % due a rate normalized plot instead
            arate = [];
            irate = [];
            for zzk = 1:size(zz1,1)
              ma = [];
              for uk = 1:2  % do over both conditions so we can keep on same scale
               ma = [ma max(RATE_att{uk}(zz1(zzk),:))];
               ma = [ma max(RATE_ign{kk}(zz1(zzk),:))];
              end
              maxo = max(ma);
              arate = [arate ; (RATE_att{kk}(zz1(zzk),:)/maxo)];
              irate = [irate ; (RATE_ign{kk}(zz1(zzk),:)/maxo)];
            end
            UU{nn,kk}{1} = nanmean(arate);
            SU{nn,kk}{1} = nanstd(arate)/sqrt(size(arate,1));
            UU{nn,kk}{2} = nanmean(irate);
            SU{nn,kk}{2} = nanstd(irate)/sqrt(size(irate,1));
          else
            UU{nn,kk}{1} = nanmean(RATE_att{kk}(zz1,:));
            SU{nn,kk}{1} = nanstd(RATE_att{kk}(zz1,:))/sqrt(length(RATE_att{kk}(zz1,:)));
            UU{nn,kk}{2} = nanmean(RATE_ign{kk}(zz1,:));
            SU{nn,kk}{2} = nanstd(RATE_ign{kk}(zz1,:))/sqrt(length(RATE_ign{kk}(zz1,:)));
          end
      end
      %*******
      if (kk == 1)
        spikeplot.plot_attention_traces(TT{nn,kk},UU{nn,kk},SU{nn,kk},LockWin,LockInt,Event);
      else
        spikeplot.plot_attention_traces(TT{nn,kk},UU{nn,kk},SU{nn,kk},LockWin2,LockInt2,Event2);    
      end
      ylabel('Rate');
      V = axis;
      VL = 0;
      VM = V(4);
      if (normo)
          VM = 0.8;
      end
      axis([V(1) V(2) VL VM]);
      aa = [LockInt2(1),LockInt2(1),LockInt2(2),LockInt2(2)];
      aa = aa + 0.001;
      bb = [VL,VM,VM,VL];
      fill(aa,bb,[0.4,0.4,0.7],'FaceAlpha',0.3,'Linestyle','none');
      if (1) %(kk == 1)
        V = axis;
        if (nn == 1)
           h4 = text(V(1)+0.1*(V(2)-V(1)),0.9*VM,sprintf('MTC'));
        end
        if (nn == 2)
           h4 = text(V(1)+0.1*(V(2)-V(1)),0.9*VM,sprintf('MT'));
        end
        if (nn == 3)
           h4 = text(V(1)+0.1*(V(2)-V(1)),0.9*VM,sprintf('Pooled'));
        end
        set(h4,'Fontsize',12);
        h4 = text(V(1)+0.1*(V(2)-V(1)),0.8*VM,sprintf('N=%d',length(zz1)));
        set(h4,'Fontsize',12);
      end
      if (kk == 1)
          ylabel('SacOn Rate');
      else
          ylabel('StimOn Rate');
      end
    end
  end
end
%**************
Mcolo = [0,0.7,0.7];
Ecolo = [0.25,0.5,0.25];
 
if (1) || (showall == 1)   %****** now plot Sensitivity Scatter plot ******
 Mcolo = [0,0.7,0.7];
 Ecolo = [0.25,0.5,0.25];
 for zk = 1:1  % SacOn vs StimOn  
   Anim = VONBASE{zk}(:,1); % MTC, MT, or Laminar?
   if (zk == 1)
       LokInt = LockInt;
   else
       LokInt = LockInt2;
   end
   hf = figure(22+zk);
   set(hf,'Position',[(200+(zk*50)) 150 1000 800]);
   FntSize = 10;
   %*************
   for nn = 1:3
    zz1 = find( Anim == nn);  % MTC 1, MT Tung 2, Pooled
    for kk = 1:3 
      ha = subplot('Position',[(0.1+(0.3*(nn-1))) (0.1+(0.3*(kk-1))) 0.20 0.20]);
      if (nn == 3)
        %****** form histograms of effects via AI computations
        zz1 = find( Anim == 1); % MTC
        zz2 = find( Anim == 2);
        switch kk
            case 1
              topname = 'Stats for FR between animals';  
              ai1 = (FR{zk}(zz1,3) - FR{zk}(zz1,5) ) ./ (FR{zk}(zz1,3) + FR{zk}(zz1,5));
              ai2 = (FR{zk}(zz2,3) - FR{zk}(zz2,5) ) ./ (FR{zk}(zz2,3) + FR{zk}(zz2,5));                
            case 2
              topname = 'Stats for MIC between animals';  
              ai1 = (MIC{zk}(zz1,3) - MIC{zk}(zz1,5) ) ./ (MIC{zk}(zz1,3) + MIC{zk}(zz1,5));
              ai2 = (MIC{zk}(zz2,3) - MIC{zk}(zz2,5) ) ./ (MIC{zk}(zz2,3) + MIC{zk}(zz2,5));                            
            case 3
              topname = 'Stats for FF between animals';
              if (0)
                ai1 = (FS{zk}(zz1,3) - FS{zk}(zz1,5) ) ./ (FS{zk}(zz1,3) + FS{zk}(zz1,5));
                ai2 = (FS{zk}(zz2,3) - FS{zk}(zz2,5) ) ./ (FS{zk}(zz2,3) + FS{zk}(zz2,5));                
              else
                ai1 = (FF{zk}(zz1,3) - FF{zk}(zz1,5) ) ./ (FF{zk}(zz1,3) + FF{zk}(zz1,5));
                ai2 = (FF{zk}(zz2,3) - FF{zk}(zz2,5) ) ./ (FF{zk}(zz2,3) + FF{zk}(zz2,5));
              end
        end
        %********
        [amed1,amed2] = compute_ai_print(ai1,ai2,topname);  % print results
        %********
        aitot = [ai1 ; ai2];
        vai = nanstd(aitot);
        if (kk == 3)
          minvai = -floor(300*vai)/100;
          maxvai = -minvai;
          step = (maxvai-minvai)/10;
        else
          minvai = -floor(300*vai)/100;
          maxvai = -minvai;
          step = (maxvai-minvai)/16;      
        end
        vx = minvai:step:maxvai;
        vy1 = hist(ai1,vx);
        vy1 = vy1 / sum(vy1);
        vy2 = hist(ai2,vx);
        vy2 = vy2 / sum(vy2);
        %***
        maxo = max(max(vy1),max(vy2));
        plot(vx,vy1,'k-','LineWidth',2,'Color',Ecolo); hold on;
        plot(vx,vy2,'k-','LineWidth',2,'Color',Mcolo); hold on;
        plot([0,0],[0,(1.2*maxo)],'k-');
        plot([amed1,amed1],[1.05*maxo,1.2*maxo],'k-','Color',Ecolo,'LineWidth',2);
        plot([amed2,amed2],[1.05*maxo,1.2*maxo],'k-','Color',Mcolo,'LineWidth',2);
        xlabel('AI index');
        ylabel('Prob');
        axis([minvai maxvai 0 (maxo*1.2)]);
        %**********
      else
        %****
        switch kk
          case 1  
            WithStd = 0; WithLog = 1;
            if (1)
              plot_von_scatter(FR{zk}(zz1,:),LokInt,WithLog,WithStd,[0.1,300],1);
              h1 = xlabel('Away FR','FontSize',FntSize);
              h2 = ylabel('Towards FR','FontSize',FntSize);
            else
              plot_von_scatter(AUC{zk}(zz1,:),LokInt,WithLog,WithStd,[0.3,1],1);
              h1 = xlabel('Away AUC','FontSize',FntSize);
              h2 = ylabel('Towards AUC','FontSize',FntSize);
            end
          case 2
            WithStd = 0; WithLog = 0;
            if (1)
              plot_von_scatter(MIC{zk}(zz1,:),LokInt,WithLog,WithStd,[0.0,1],1);
              h1 = xlabel('Away MIC','FontSize',FntSize);
              h2 = ylabel('Towards MIC','FontSize',FntSize);
            else
              plot_von_scatter(MIC2{zk}(zz1,:),LokInt,WithLog,WithStd,[0.0,1],1);
              h1 = xlabel('Away MIC','FontSize',FntSize);
              h2 = ylabel('Towards MIC','FontSize',FntSize);     
            end
          case 3
            WithStd = 0; WithLog = 0;
            if (0)
              plot_von_scatter(FS{zk}(zz1,:),LokInt,WithLog,WithStd,[0.0,2],1);
              h1 = xlabel('Away FS','FontSize',FntSize);
              h2 = ylabel('Towards FS','FontSize',FntSize);                
            else
              plot_von_scatter(FF{zk}(zz1,:),LokInt,WithLog,WithStd,[0.0,2],1);
              h1 = xlabel('Away FF','FontSize',FntSize);
              h2 = ylabel('Towards FF','FontSize',FntSize);
            end
      end
      set(h1,'Color',[0,0,1]);
      set(h2,'Color',[1,0,0]);
      if (kk == 3) 
        V = axis;
        va = V(1)+0.1*(V(2)-V(1));
        vb = V(3)+0.9*(V(4)-V(3));    
        vb2 = V(3)+0.8*(V(4)-V(3));
        if (nn == 1)
           h4 = text(va,vb,sprintf('MTC'));
        end
        if (nn == 2)
           h4 = text(va,vb,sprintf('MT'));
        end
        if (nn == 3)
           h4 = text(va,vb,sprintf('MTLam'));
        end
        set(h4,'Fontsize',FntSize);
        h4 = text(va,vb2,sprintf('N=%d',length(zz1)));
        set(h4,'Fontsize',FntSize);
      end
      set(gca,'FontSize',FntSize);
      %********
    end
    end
   end
  end
end


if (0) || (showall == 1)  % plotting the mean change in VonMises Fit
  yb = 0.85/(length(VONFIT)-1);
  cmap = [[0.80,0.80,0.40]; ...
        [0.87,0.70,0.43]; ...
        [0.93,0.61,0.46]; ...
        [1.00,0.50,0.50]];
  %***** load in Von Mises parameters ******
  acurve = cell(3,2);
  acurve1 = cell(3,2);
  acurve2 = cell(3,2);
  ucurve = cell(3,2);
  ucurve1 = cell(3,2);
  ucurve2 = cell(3,2);
  normo = 1;
  %********
  Anim = VONFIT{1}(:,1); % MTC, MT, or Laminar? 
  hf = figure(24);
  set(hf,'Position',[400 150 1200 800]);
  %*************
  for nn = 1:3
    zz1 = find( (Anim == nn) & (R2List == 1) );  % MTC 1, MT Tung 2, MT Laminar 3
    for kk = 1:length(VONFIT)
      ha = subplot('Position',[(0.1+(0.3*(nn-1))) (0.1+(0.5*(kk-1))) 0.25 0.35]); 
      if nn == 3
        acurve{3,kk} = 0.5*(acurve{1,kk}+acurve{2,kk});
        acurve1{3,kk} = 0.5*(acurve1{1,kk}+acurve1{2,kk});
        acurve2{3,kk} = 0.5*(acurve2{1,kk}+acurve2{2,kk});
        ucurve{3,kk} = 0.5*(ucurve{1,kk}+ucurve{2,kk});
        ucurve1{3,kk} = 0.5*(ucurve1{1,kk}+ucurve1{2,kk});
        ucurve2{3,kk} = 0.5*(ucurve2{1,kk}+ucurve2{2,kk});
        zz1 = [1,2];
      else
        if (kk == 1) % length(VONFIT))
          xlab = 1;
        else
          xlab = 0;
        end
        Ibo = VONFIT{kk}(zz1,:);
        %******* average the raw curves, but peak normalize them
        if (1)
          acurves = [];
          ucurves = [];
          maxos = [];
          x = (0:360)*(pi/180);
          for ck = 1:size(Ibo,1) 
            amu = Ibo(ck,3:6);
            bmu = Ibo(ck,7:10);
            if (1)
              amu(4) = pi;  % fix to standard center
              bmu(4) = pi;  % or could take the average as center...
            else  
              dpref = amu(4)-bmu(4);
              if (dpref > pi)
                dpref = (2*pi)-dpref;
              end
              if (dpref < -pi)
                dpref = dpref + (2*pi);
              end
              amu(4) = pi + (dpref/2);
              bmu(4) = pi - (dpref/2);
            end
            acur = VonMises.vonmises(x,amu); 
            bcur = VonMises.vonmises(x,bmu);
            if (normo)
              maxo = 0.5*( max(acur) + max(bcur) );
              % maxo = min(nanmax(acur),nanmax(bcur));
              acur = acur / maxo;  % norm to peak
              bcur = bcur / maxo;  % norm to peak
            end
            acurves = [acurves ; acur];
            ucurves = [ucurves ; bcur];
            maxos = [maxos ; maxo];
          end
          if (normo)
              mmaxo = 1;
          else
              mmaxo = nanmean(maxos);
          end
          ucurve{nn,kk} = nanmean(ucurves) * mmaxo; 
          suu = (2 * nanstd(ucurves) * mmaxo)/sqrt(length(zz1));
          ucurve1{nn,kk} = ucurve{nn,kk} + suu;
          ucurve2{nn,kk} = ucurve{nn,kk} - suu;
          acurve{nn,kk} = nanmean(acurves) * mmaxo;
          suu = (2 * nanstd(acurves) * mmaxo)/sqrt(length(zz1));
          acurve1{nn,kk} = acurve{nn,kk} + suu;
          acurve2{nn,kk} = acurve{nn,kk} - suu;
          %******
        else
          %*** average in the von mises parameter space instead
          %****** this will not work for k == 3 anymore
          for kz = 1:4
            uu(kz) = nanmedian(Ibo(:,6+kz));
            suu(kz) = nanstd(Ibo(:,(6+kz)))/sqrt(length(Ibo(:,(6+kz))));
            if (kz < 3)
              aimod = nanmedian( (Ibo(:,2+kz) - Ibo(:,6+kz)) ./ (Ibo(:,2+kz)+Ibo(:,6+kz)) );
              aa(kz) = uu(kz)*((1+aimod)/(1-aimod));
              saa(kz) = nanstd(Ibo(:,(2+kz)))/sqrt(length(Ibo(:,(2+kz))));
            else
              aimod = nanmedian( (Ibo(:,2+kz) - Ibo(:,6+kz)) );  % fit as mean change
              aa(kz) = uu(kz) + aimod;
              saa(kz) = nanstd(Ibo(:,(2+kz)))/sqrt(length(Ibo(:,(2+kz))));    
            end
          end
          %****
          x = (0:360)*(pi/180);
          uu(4) = pi;
          aa(4) = pi;
          suu(4) = 0;
          saa(4) = 0;
          ucurve{nn,kk} = VonMises.vonmises(x,uu);
          ucurve1{nn,kk} = VonMises.vonmises(x,uu+suu);
          ucurve2{nn,kk} = VonMises.vonmises(x,uu-suu);
          %*****
          acurve{nn,kk} = VonMises.vonmises(x,aa);
          acurve1{nn,kk} = VonMises.vonmises(x,aa+saa);
          acurve2{nn,kk} = VonMises.vonmises(x,aa-saa);
        end
      end
      %**********
      fcolor = cmap(kk,:);
      fcolor2 = fcolor;
      fcolor2(2:3) = fcolor(2:3)*0.5;
      %******
      xx = [x fliplr(x)];
      yy = [ucurve1{nn,kk} fliplr(ucurve2{nn,kk})];
      fill(xx,yy,[0.5,0.5,1.0],'FaceAlpha',0.3,'Linestyle','none'); hold on;  
      h = plot(x,ucurve{nn,kk},'b.'); hold on;
      set(h,'LineWidth',2);
      %******
      xx = [x fliplr(x)];
      yy = [acurve1{nn,kk} fliplr(acurve2{nn,kk})];
      fill(xx,yy,fcolor,'FaceAlpha',0.3,'Linestyle','none');
      h2 = plot(x,acurve{nn,kk},'r-');
      set(h2,'LineWidth',2);
      set(h2,'Color',fcolor2);
      axis tight;
      if ~normo
        plot([0,(2*pi)],[20,20],'k--');
        hh = ylabel('Rate (sp/s)');
      else
        hh = ylabel('Norm Rate (sp/s)');     
      end
      set(hh,'Fontsize',16);
      if (xlab)
        hh = xlabel('Motion Direction');
        set(hh,'Fontsize',16);
      end
      if (1) %(kk == 1)
        V = axis;
        axis([V(1) V(2) 0 V(4)]);
        if (nn == 1)
           h4 = text((1.5*pi),0.9*V(4),sprintf('MTC'));
        end
        if (nn == 2)
           h4 = text((1.5*pi),0.9*V(4),sprintf('MT'));
        end
        if (nn == 3)
           h4 = text((1.5*pi),0.9*V(4),sprintf('Pooled'));
        end
        set(h4,'Fontsize',12);
        h4 = text((1.5*pi),0.8*V(4),sprintf('N=%d',length(zz1)));
        set(h4,'Fontsize',12);
      end
      if (kk == 1)
          ylabel('SacOn Rate');
      else
          ylabel('StimOn Rate');
      end
    end
  end
end

if (0) || (showall == 1)    %****** now plot VonMises Params Scatter plot ******
  for zk = 1:1  % SacOn vs StimOn  
   Anim = VONBASE{zk}(:,1); % MTC, MT, or Laminar?
   if (zk == 1)
       LokInt = LockInt;
   else
       LokInt = LockInt2;
   end
   hf = figure(25+zk);
   set(hf,'Position',[(500+(zk*50)) 150 1000 800]);
   FntSize = 10;
   %*************
   %*************
   for nn = 1:3
    zz1 = find( (Anim == nn) & (R2List == 1) );  % MTC 1, MT Tung 2, Pooled
    for kk = 1:3 
      ha = subplot('Position',[(0.1+(0.3*(nn-1))) (0.1+(0.3*(kk-1))) 0.20 0.20]);
      if (nn == 3)
        %****** form histograms of effects via AI computations
        zz1 = find( (Anim == 1) & (R2List == 1) ); % MTC
        zz2 = find( (Anim == 2) & (R2List == 1) );
        switch kk
            case 1
              topname = 'Stats for BASE between animals';  
              ai1 = (VONBASE{zk}(zz1,3) - VONBASE{zk}(zz1,5) ) ./ ...
                    (VONBASE{zk}(zz1,3) + VONBASE{zk}(zz1,5));
              ai2 = (VONBASE{zk}(zz2,3) - VONBASE{zk}(zz2,5) ) ./ ...
                    (VONBASE{zk}(zz2,3) + VONBASE{zk}(zz2,5));                
            case 2
              topname = 'Stats for GAIN between animals';  
              ai1 = (VONGAIN{zk}(zz1,3) - VONGAIN{zk}(zz1,5) ) ./ ...
                    (VONGAIN{zk}(zz1,3) + VONGAIN{zk}(zz1,5));
              ai2 = (VONGAIN{zk}(zz2,3) - VONGAIN{zk}(zz2,5) ) ./ ...
                    (VONGAIN{zk}(zz2,3) + VONGAIN{zk}(zz2,5));                                            
            case 3
              topname = 'Stats for WIDTH between animals'; 
              ai1 = (VONKAP{zk}(zz1,3) - VONKAP{zk}(zz1,5) ) ./ ...
                    (VONKAP{zk}(zz1,3) + VONKAP{zk}(zz1,5));
              ai2 = (VONKAP{zk}(zz2,3) - VONKAP{zk}(zz2,5) ) ./ ...
                    (VONKAP{zk}(zz2,3) + VONKAP{zk}(zz2,5));                                            
        end
        %********
        [amed1,amed2] = compute_ai_print(ai1,ai2,topname);  % print results
        %********** plot histograms
        aitot = [ai1 ; ai2];
        vai = nanstd(aitot);
        minvai = -floor(300*vai)/100;
        maxvai = -minvai;
        step = (maxvai-minvai)/12;
        vx = minvai:step:maxvai;
        vy1 = hist(ai1,vx);
        vy1 = vy1 / sum(vy1);
        vy2 = hist(ai2,vx);
        vy2 = vy2 / sum(vy2);
        %***
        maxo = max(max(vy1),max(vy2));
        plot(vx,vy1,'k-','LineWidth',2,'Color',Ecolo); hold on;
        plot(vx,vy2,'k-','LineWidth',2,'Color',Mcolo); hold on;
        plot([0,0],[0,(1.2*maxo)],'k-');
        plot([amed1,amed1],[1.05*maxo,1.2*maxo],'k-','Color',Ecolo,'LineWidth',2);
        plot([amed2,amed2],[1.05*maxo,1.2*maxo],'k-','Color',Mcolo,'LineWidth',2);
        xlabel('AI index');
        ylabel('Prob');
        axis([minvai maxvai 0 (maxo*1.2)]);
        %**********
      else
        %****
        switch kk
          case 1  
            WithStd = 0; WithLog = 1; Range = [0.1,200];
            plot_von_scatter(VONBASE{zk}(zz1,:),LokInt,WithLog,WithStd,Range,1);
            h1 = xlabel('Away Base','FontSize',FntSize);
            h2 = ylabel('Towards Base','FontSize',FntSize);
          case 2
            WithStd = 0; WithLog = 1; Range = [5,200];
            plot_von_scatter(VONGAIN{zk}(zz1,:),LokInt,WithLog,WithStd,Range,1);
            h1 = xlabel('Away Gain','FontSize',FntSize);
            h2 = ylabel('Towards Gain','FontSize',FntSize);
          case 3
            WithStd = 0; WithLog = 0; Range = [0,180]; % Range = [0.01,10];
            plot_von_scatter(VONKAP{zk}(zz1,:),LokInt,WithLog,WithStd,Range,1);
            h1 = xlabel('Away Kappa','FontSize',FntSize);
            h2 = ylabel('Towards Kappa','FontSize',FntSize);
        end
        set(h1,'Color',[0,0,1]);
        set(h2,'Color',[1,0,0]);
      
        if (kk == 3) 
          if (WithLog)
            Range = log10(Range);
          end
          va = Range(1)+0.1*(Range(2)-Range(1));
          vb = Range(1)+0.9*(Range(2)-Range(1));    
          vb2 = Range(1)+0.8*(Range(2)-Range(1));
          if (WithLog)
            va = 10^va;
            vb = 10^vb;
            vb2 = 10^vb2;
          end
          if (nn == 1)
           h4 = text(va,vb,sprintf('MTC'));
          end
          if (nn == 2)
           h4 = text(va,vb,sprintf('MT'));
          end
          if (nn == 3)
           h4 = text(va,vb,sprintf('MTLam'));
          end
          set(h4,'Fontsize',FntSize);
          h4 = text(va,vb2,sprintf('N=%d',length(zz1)));
          set(h4,'Fontsize',FntSize);
        end
        set(gca,'FontSize',FntSize);
        %********
      end
    end
   end
  end
end

if (0) || (showall == 1)   %****** now non-parametric tuning fits ******
  normo = 1;
  Anim = VONBASE{zk}(:,1); % MTC, MT, or Laminar?  
  hf = figure(32);
  set(hf,'Position',[100 150 1200 800]);
  TT = cell(3,2);  UU = cell(3,2);  SU = cell(3,2); % 1-3 animal, 1-2 timeint
  AIBASE = cell(3,2); AIAMP = cell(3,2);  AIWIDTH = cell(3,2);
  %*******
  LabelA = 'SacOn Rate(All)';
  LabelB = 'FeatG Rate(All)';
  %*******
  for nn = 1:3  
    for kk = 1:2 
      if (kk == 1)
          zz1 = find( (Anim == nn) );  % MTC 1, MT Tung 2, MT Laminar 3
          HOLD = NPFIT;
          
      else
          zz1 = find( (Anim == nn) );  % MTC 1, MT Tung 2, MT Laminar 3
          HOLD = NPFEAT;
      end
      ha = subplot('Position',[(0.1+(0.3*(nn-1))) (0.1+(0.5*(kk-1))) 0.20 0.35]);
      if nn == 3
          %****** pool the results from the previous two and plot them
          TT{3,kk}{1} = 0.5*(TT{1,kk}{1} + TT{2,kk}{1});
          TT{3,kk}{2} = 0.5*(TT{1,kk}{2} + TT{2,kk}{2});
          UU{3,kk}{1} = 0.5*(UU{1,kk}{1} + UU{2,kk}{1});
          SU{3,kk}{1} = 0.5*(SU{1,kk}{1} + SU{2,kk}{1});
          UU{3,kk}{2} = 0.5*(UU{1,kk}{2} + UU{2,kk}{2});
          SU{3,kk}{2} = 0.5*(SU{1,kk}{2} + SU{2,kk}{2});
          AIBASE{3,kk} = [AIBASE{1,kk} ; AIBASE{2,kk}];
          AIAMP{3,kk} = [AIAMP{1,kk} ; AIAMP{2,kk}];
          zz1 = [1,2];  % two animals
          %******        
      else
          LockO = [1,size(HOLD{1},2)];
          TT{nn,kk}{1} = 1:size(HOLD{1},2);
          TT{nn,kk}{2} = TT{nn,kk}{1};
          if (normo == 1) % due a rate normalized plot instead
               arate = [];
               irate = [];
               aibase = [];
               aiamp = [];
               aiwid = [];
               for zzk = 1:size(zz1,1)
                  ma = [];
                  mb = [];
                  for uk = 1:2  % do over both conditions so we can keep on same scale
                     ma = [ma max(HOLD{uk}(zz1(zzk),:))];
                     mb = [mb min(HOLD{uk}(zz1(zzk),:))];
                  end
                  maxo = mean(ma);
                  mino = mean(mb);
                  aa = ( HOLD{1}(zz1(zzk),:)/maxo);
                  uu = ( HOLD{2}(zz1(zzk),:)/maxo);
                  basea = nanmean(aa(1:4));
                  baseu = nanmean(uu(1:4));
                  ampa = nanmean(aa(13:16))-basea;
                  ampu = nanmean(uu(13:16))-baseu;
                  %**** non-parametric width
                  zaa = max((aa-basea),0);
                  zaa = min((zaa/ampa),1);  % bound from 0 to 1
                  zuu = max((uu-baseu),0);
                  zuu = min((zuu/ampu),1);  % bound from 0 to 1
                  awid = sum(zaa)/length(zaa);
                  uwid = sum(zuu)/length(zuu);
                  %**********
                  arate = [arate ; aa];
                  irate = [irate ; uu];
                  %*******
                  aibase = [aibase ; ((basea-baseu)/(basea+baseu))];
                  aiamp  = [aiamp  ; ((ampa-ampu)/(ampa+ampu))];
                  aiwid  = [aiwid ; ((awid-uwid)/(awid+uwid))];
               end
               UU{nn,kk}{1} = nanmean(arate);
               SU{nn,kk}{1} = nanstd(arate)/sqrt(size(arate,1));
               UU{nn,kk}{2} = nanmean(irate);
               SU{nn,kk}{2} = nanstd(irate)/sqrt(size(irate,1));
               AIBASE{nn,kk} = aibase;
               AIAMP{nn,kk} = aiamp;
               AIWID{nn,kk} = aiwid;
          else
               UU{nn,kk}{1} = nanmean(HOLD{1}(zz1,:));
               SU{nn,kk}{1} = nanstd(HOLD{1}(zz1,:))/sqrt(size(HOLD{1}(zz1,:),1));
               UU{nn,kk}{2} = nanmean(HOLD{2}(zz1,:));
               SU{nn,kk}{2} = nanstd(HOLD{2}(zz1,:))/sqrt(size(HOLD{2}(zz1,:),1));
          end  
      end    
      %*******
      spikeplot.plot_attention_traces(TT{nn,kk},UU{nn,kk},SU{nn,kk},LockO,[],[]);
      ylabel('Rate');
      V = axis;
      VL = 0;
      VM = V(4);
      if (normo)
          VM = 1.0;
      end
      axis([V(1) V(2) VL VM]);
      if (1) %(kk == 1)
        V = axis;
        if (nn == 1)
           h4 = text(V(1)+0.1*(V(2)-V(1)),0.9*VM,sprintf('MTC'));
        end
        if (nn == 2)
           h4 = text(V(1)+0.1*(V(2)-V(1)),0.9*VM,sprintf('MT'));
        end
        if (nn == 3)
           h4 = text(V(1)+0.1*(V(2)-V(1)),0.9*VM,sprintf('Pooled'));
        end
        set(h4,'Fontsize',12);
        h4 = text(V(1)+0.1*(V(2)-V(1)),0.8*VM,sprintf('N=%d',length(zz1)));
        set(h4,'Fontsize',12);
      end
      if (kk == 1)
          ylabel(LabelA); 
      else
          ylabel(LabelB); 
      end
    end
  end
end
%**** report some stats then *******
if (0) || (showall == 1)
  disp(' ');
  disp('Non-parametric stats (all units) ....');
  for zk = 1:2
    disp('  ');
    if (zk == 1)
      disp('... for towards vs away conditions');
    else
      disp('... for feature gain, all way, pref or non-pref');
    end
%     %******* report non-parametric results0
    compute_ai_print(AIBASE{1,zk},AIBASE{2,zk},'Baseline rates ');
    compute_ai_print(AIAMP{1,zk},AIAMP{2,zk},'Amplitudes ');
    compute_ai_print(AIWID{1,zk},AIWID{2,zk},'Widths ');
    %*******
  end
  %*****
  if (0)
      disp(' ');
      disp('Plotting Width Scatter');
      if (0) %substitute in Von Mises Widths instead
         zz = find( R2List == 1);
         NPWID = [VONKAP{1}(zz,1) VONKAP{1}(zz,3) 0.5*VONKAP{1}(zz,4) VONKAP{1}(zz,5) 0.5*VONKAP{1}(zz,6)];
         zz = 1:length(NPWID);
         wmino = 0;
         wmaxo = 180;
      else
         zz = find( (NPWID(:,4) > 0) & (NPWID(:,2)> 0) & (NPWID(:,4)<1) & (NPWID(:,2)<1) );
         wmino = 0;
         wmaxo = 1;
      end
      zz1 = find( NPWID(zz,1) == 1); % MTC
      zz2 = find( NPWID(zz,1) == 2); % MT
      comp = NPWID(zz,2) - NPWID(zz,4);
      wido = 1.0*(NPWID(zz,3)+NPWID(zz,5));  %95 perc conf for Gaussian (conservative)
      abo = abs(comp) ./ wido;
      sigo = (abo > 1.96); 
      iz = find(sigo > 0);
      nzz1 = zz( find( (sigo == 0) & (NPWID(zz,1) == 1)) );
      szz1 = zz( find( (sigo == 1) & (NPWID(zz,1) == 1)) );
      nzz2 = zz( find( (sigo == 0) & (NPWID(zz,1) == 2)) );
      szz2 = zz( find( (sigo == 1) & (NPWID(zz,1) == 2)) );
      disp(sprintf('%d of %d gave measurable values',length(zz),length(NPWID)));
      disp(sprintf('%d of %d gave significant mod',length(iz),length(zz)));
      %**********
      hf = figure(50);
      set(hf,'Position',[600 400 1200 600]);
      %******
      subplot('Position',[0.1 0.15 0.375 0.75]);
      plot(NPWID(nzz2,4),NPWID(nzz2,2),'ko','Color',Mcolo,'Markersize',5); hold on;
      plot(NPWID(szz2,4),NPWID(szz2,2),'k.','Color',Mcolo,'Markersize',10); 
      plot(NPWID(nzz1,4),NPWID(nzz1,2),'ko','Color',Ecolo,'Markersize',5);
      plot(NPWID(szz1,4),NPWID(szz1,2),'k.','Color',Ecolo,'Markersize',10);
      axis([wmino wmaxo wmino wmaxo]);
      plot([wmino,wmaxo],[wmino,wmaxo],'k-');
      xlabel('Away NP Width','Color',[0,0,1],'Fontsize',18);
      ylabel('Towards NP Width','Color',[1,0,0],'Fontsize',18);
      %***
      vmin = -0.6;
      vmax = 0.6;
      vx = vmin:0.04:vmax;
      %*****
      subplot('Position',[0.6 0.15 0.35 0.35]);
      nai = (NPWID(nzz2,2)-NPWID(nzz2,4)) ./ (NPWID(nzz2,2)+NPWID(nzz2,4));
      sai = (NPWID(szz2,2)-NPWID(szz2,4)) ./ (NPWID(szz2,2)+NPWID(szz2,4));
      ys = hist(sai,vx);
      yn = hist(nai,vx);
      yy = ys+yn;
      maxo = max(yy);
      plot(vx,ys,'k.-','Linewidth',2); hold on;
      plot(vx,yy,'k.-','Color',Mcolo,'Linewidth',2);
      axis([vmin vmax 0 maxo]);
      plot([0,0],[0,maxo],'k-');
      xlabel('AI NP Width','Fontsize',18);
      ylabel('Number of Cells','Color',Mcolo,'Fontsize',18);
      text(0.1,0.8*maxo,[sprintf('%4.1f',100*length(sai)/...
          (length(sai)+length(nai))),'% sig'],'Color',Mcolo,'Fontsize',12);
      ai = [sai ; nai];
      medai = median(ai);
      [p,h] = signrank(ai);
      text(0.1,0.9*maxo,[sprintf('Med %5.3f(p=%5.4f)',((1+medai)/(1-medai)),p)],...
                 'Color',Mcolo,'Fontsize',12);
      %*********
      vmin = -0.6;
      vmax = 0.6;
      vx = vmin:0.06:vmax;
      subplot('Position',[0.6 0.55 0.35 0.35]);
      nai = (NPWID(nzz1,2)-NPWID(nzz1,4)) ./ (NPWID(nzz1,2)+NPWID(nzz1,4));
      sai = (NPWID(szz1,2)-NPWID(szz1,4)) ./ (NPWID(szz1,2)+NPWID(szz1,4));
      ys = hist(sai,vx);
      yn = hist(nai,vx);
      yy = ys+yn;
      maxo = max(yy);
      plot(vx,ys,'k.-','Linewidth',2); hold on;
      plot(vx,yy,'k.-','Color',Ecolo,'Linewidth',2);
      axis([vmin vmax 0 maxo]);
      plot([0,0],[0,maxo],'k-');
      % xlabel('AI NP Width','Fontsize',18);
      ylabel('Number of Cells','Color',Ecolo,'Fontsize',18);
      text(0.1,0.8*maxo,[sprintf('%4.1f',100*length(sai)/...
          (length(sai)+length(nai))),'% sig'],'Color',Ecolo,'Fontsize',12);
      ai = [sai ; nai];
      medai = median(ai);
      [p,h] = signrank(ai);
      text(0.1,0.9*maxo,[sprintf('Med %5.3f(p=%5.4f)',(1+medai)/(1-medai),p)],...
                 'Color',Ecolo,'Fontsize',12);
      %*******
  end
  %**********
end
%**************

%****** quick look at NCORR stats
if (0)
NBOUND = 0.5*(NCORR(:,5)+NCORR(:,6));
NSTEP = 0.05;
NSMO = 0.10;
%*********
hf = figure(40);
set(hf,'Position',[100 100 1200 800]);
%*** only include pairs where the corr had reasonable bounds on its estimate        
for plt = 1:4 %7
 %******* 
 if (plt == 1)
    subplot('Position',[0.1 0.6 0.35 0.35]);    
    kzz = 1:length(NCORR);
    Tlabel = 'MT Noise Correlations';
 else
    if (0) 
       layer = floor((plt-2)/2);
       waver = mod((plt-2),2);
       kzz = find( (NCORR(:,7)==(plt-1)) & (NCORR(:,8)==(plt-1)) );
       if ~isempty(kzz)
          subplot('Position',[ (0.5+(0.225*waver)) (0.7-(0.3*layer)) 0.2 0.2]);
          Tlabel = sprintf('N=%d L(%d)W(%d)',length(kzz),(layer+1),(waver+1));
       end
    else
       layer = plt-2;
       laya = (layer*2)+1;
       layb = (layer*2)+2;
       % laya = layb; % broad only
       % layb = laya; % narrow only
       kzz = find( ((NCORR(:,7)==laya) | (NCORR(:,7)==layb)) & ...
                ((NCORR(:,8)==laya) | (NCORR(:,8)==layb)) );
       if ~isempty(kzz)
          subplot('Position',[ 0.6 (0.7-(0.3*layer)) 0.3 0.2]);
          Tlabel = sprintf('N=%d L(%d)',length(kzz),(layer+1));
       end
    end
 end
 %*******
 if ~isempty(NCORR) & ~isempty(kzz)  % look at MT only (too few in MTC, pointless)
    vx = -1:NSTEP:1;
    xu = [];
    au = [];
    asu = [];
    iu = [];
    isu = [];
    uu = [];
    su = [];
    mu = [];
    zu = [];
    for it = 1:(length(vx)-1)
        ita = vx(it)-NSMO;
        itb = vx(it+1)+NSMO;
        zz = find( (NCORR(kzz,4) >= ita) & (NCORR(kzz,4) <= itb) ); 
        zu = [zu length(zz)];
        if (length(zz) > 1)
          di = NCORR(kzz(zz),2);
          di = di - NCORR(kzz(zz),3);  
          u = nanmean(di);
          s = nanstd(di)/sqrt(length(zz));
          m = 0.5*(nanmean(NCORR(kzz(zz),2)+NCORR(kzz(zz),3)));
          xu = [xu (0.5*(ita+itb))];
          uu = [uu u];
          su = [su s];
          mu = [mu m];
          %*****
          au = [au nanmean(NCORR(kzz(zz),2))];
          iu = [iu nanmean(NCORR(kzz(zz),3))];
          asu = [asu (nanstd(NCORR(kzz(zz),2))/sqrt(length(zz)))];
          isu = [isu (nanstd(NCORR(kzz(zz),3))/sqrt(length(zz)))];
        end
    end
    %**********
    xx = [xu fliplr(xu)];
    yy = [(au+(2*asu)) fliplr(au-(2*asu))];
    fill(xx,yy,[1,0,0],'FaceAlpha',0.2,'Linestyle','none'); hold on;
    plot(xu,au,'r-','Linewidth',2);
    yy = [(iu+(2*isu)) fliplr(iu-(2*isu))];
    fill(xx,yy,[0,0,1],'FaceAlpha',0.2,'Linestyle','none');
    plot(xu,iu,'b-','Linewidth',2);
    xlabel('Signal');
    ylabel('Corr');  
    if (plt == 1)
       axis([-1 1 -0.02 0.12]);
    else
       axis([-1 1 -0.05 0.15]);    
    end
    plot([-1,1],[0,0],'k:');
    plot([0,0],[-0.02,0.12],'k:');
    title(Tlabel);
 end
 if (plt == 1)
    subplot('Position',[0.1 0.1 0.35 0.35]);
    xx = [xu fliplr(xu)];
    yy = [(uu+(2*su)) fliplr(uu-(2*su))];
    fill(xx,yy,[0,0,0],'FaceAlpha',0.2,'Linestyle','none'); hold on;
    plot(xu,uu,'k-','Linewidth',2);
    xlabel('Signal');
    ylabel('Delta Corr.');
    axis([-1 1 -0.04 0.04]);
    plot([-1,1],[0,0],'k:');
    plot([0,0],[-0.04,0.04],'k:');
    title(sprintf('N = %d pairs',length(NBOUND)));
 end
end
end
%%
%********* Laminar Depth Scatter Plots
if (1)
USENP = 0;  % use not parametric (all units) or use Von Mises (R2>thresh)
MINRATE = 1;
wthreshmin = Wthreshmin;
wthreshmax = Wthreshmax;
wthresh = mean([wthreshmin,wthreshmax]);
if (1) || (showall == 1)
    depth = CSD-NDepth;
    if (USENP) % use all units with non-parametric
       zz = find( ~isnan(depth) & (NPMIN > MINRATE) );  % all units with laminar info        
       ndepth = depth(zz);
       nduration = Duration(zz); 
       narrow = find( nduration < wthreshmin);
       broad = find( nduration > wthreshmax);
       %****** non-parametric forms
       adelta = NPAMP(zz,3) - NPAMP(zz,5);  % gain change
       aigain = adelta ./ (NPAMP(zz,3) + NPAMP(zz,5));  %norm AI
       bdelta = NPBASE(zz,3) - NPBASE(zz,5);  % base change
       aibase = bdelta ./ (NPBASE(zz,3) + NPBASE(zz,5));  %norm AI
       wdelta = NPWID(zz,3) - NPWID(zz,5);  % width change
       aiwid = wdelta ./ (NPWID(zz,3) + NPWID(zz,5));  %norm AI 
       %****** non-parametric forms
       adelta = MIC{1}(zz,3) - MIC{1}(zz,5);  % gain change
       aimic = adelta ./ (MIC{1}(zz,3) + MIC{1}(zz,5));  %norm AI
       bdelta = FF{1}(zz,3) - FF{1}(zz,5);  % base change
       aiff = bdelta ./ (FF{1}(zz,3) + FF{1}(zz,5));  %norm AI
       cdelta = FR{1}(zz,3) - FR{1}(zz,5);  % base change
       airate = cdelta ./ (FR{1}(zz,3) + FR{1}(zz,5));  %norm AI
    else  % use well fit Von Mises
       zz = find( ~isnan(depth) & (R2List == 1) & (NPMIN > MINRATE));
       ndepth = depth(zz);
       nduration = Duration(zz);
       narrow = find( nduration < wthreshmin);
       broad = find( nduration > wthreshmax);
       %**** von mises fits
       adelta = VONGAIN{1}(zz,3) - VONGAIN{1}(zz,5);  % gain change
       aigain = adelta ./ (VONGAIN{1}(zz,3) + VONGAIN{1}(zz,5));  %norm AI
       bdelta = VONBASE{1}(zz,3) - VONBASE{1}(zz,5);  % base change
       aibase = bdelta ./ (VONBASE{1}(zz,3) + VONBASE{1}(zz,5));  %norm AI
       wdelta = VONKAP{1}(zz,3) - VONKAP{1}(zz,5);  % width change
       aiwid = wdelta ./ (VONKAP{1}(zz,3) + VONKAP{1}(zz,5));  %norm AI
       adelta = MIC{1}(zz,3) - MIC{1}(zz,5);  % gain change
       aimic = adelta ./ (MIC{1}(zz,3) + MIC{1}(zz,5));  %norm AI
       bdelta = FF{1}(zz,3) - FF{1}(zz,5);  % base change
       aiff = bdelta ./ (FF{1}(zz,3) + FF{1}(zz,5));  %norm AI
       cdelta = FR{1}(zz,3) - FR{1}(zz,5);  % base change
       airate = cdelta ./ (FR{1}(zz,3) + FR{1}(zz,5));  %norm AI
    end
    
    %******* plot parameters
    msize = 7;
    wmsize = 10;
    ncolo = [0.6,0.6,0];
    bcolo = [0,0.6,0];
    aivx = -0.56:0.08:0.56;
    wvx = 0:1:ceil(2*wthresh);
    %**************
    
    if (1) % plot tuning function parameters
        %********
        hf = figure;
        set(hf,'Position',[400 200 1200 800]);
        %******
        znduration = nduration + (rand(size(nduration))-0.5);
        scat_panel(0.1,ndepth,znduration,narrow,broad,ncolo,bcolo,wmsize,wvx,0);
        xlabel('Duration');    
        %******
        scat_panel(0.325,ndepth,aigain,narrow,broad,ncolo,bcolo,msize,aivx,1);
        xlabel('Gain');
        %******
        scat_panel(0.55,ndepth,aibase,narrow,broad,ncolo,bcolo,msize,aivx,1);
        xlabel('Base');
        %******
        scat_panel(0.775,ndepth,aiwid,narrow,broad,ncolo,bcolo,msize,aivx,1);
        xlabel('Width');
    end
    if (1) % plot tuning function parameters
        %********
        hf = figure;
        set(hf,'Position',[500 200 1200 800]);
        %******
        znduration = nduration + (rand(size(nduration))-0.5);
        scat_panel(0.1,ndepth,znduration,narrow,broad,ncolo,bcolo,wmsize,wvx,0);
        xlabel('Duration');    
        %******
        scat_panel(0.325,ndepth,aimic,narrow,broad,ncolo,bcolo,msize,aivx,1);
        xlabel('MIC');
        %******
        scat_panel(0.55,ndepth,aiff,narrow,broad,ncolo,bcolo,msize,aivx,1);
        xlabel('FF');
        %******
        scat_panel(0.775,ndepth,airate,narrow,broad,ncolo,bcolo,msize,aivx,1);
        xlabel('Rate');
        %*****
    end
end
end
    

return;

function scat_panel(xx,ndepth,nduration,narrow,broad,ncolo,bcolo,msize,vx,correlated)
        subplot('Position',[xx 0.35 0.20 0.60]);
        plot(nduration,ndepth,'k.','Markersize',(msize/2)); hold on;
        plot(nduration(narrow),ndepth(narrow),'k.','Markersize',msize,'Color',ncolo); hold on;
        plot(nduration(broad),ndepth(broad),'k.','Markersize',msize,'Color',bcolo); hold on;
        %******* compute rolling mean over depth
        if (correlated)
            dmin = min(ndepth);
            dmax = max(ndepth);
            dbin = 100;
            minum = 10;
            nuu = [];
            nsu = [];
            buu = [];
            bsu = [];
            auu = [];
            asu = [];
            xuu = [];
            for dx = dmin:((dmax-dmin)/100):dmax
                dxa = dx-dbin;
                dxb = dx+dbin;
                xuu = [xuu ; dx];
                %******
                zz = find( (ndepth(narrow) > dxa) & (ndepth(narrow) < dxb));
                if length(zz) < minum
                    nuu = [nuu ; NaN];
                    nsu = [nsu ; NaN];
                else
                    nuu = [nuu ; nanmean(nduration(narrow(zz)))];
                    nsu = [nsu ; (nanstd(nduration(narrow(zz)))/sqrt(length(zz)))];
                end
                %*******
                zz = find( (ndepth(broad) > dxa) & (ndepth(broad) < dxb));
                if length(zz) < minum
                    buu = [buu ; NaN];
                    bsu = [bsu ; NaN];
                else
                    buu = [buu ; nanmean(nduration(broad(zz)))];
                    bsu = [bsu ; (nanstd(nduration(broad(zz)))/sqrt(length(zz)))];
                end
                %*******
                zz = find( (ndepth > dxa) & (ndepth < dxb));
                if length(zz) < minum
                    auu = [auu ; NaN];
                    asu = [asu ; NaN];
                else
                    auu = [auu ; nanmean(nduration(zz))];
                    asu = [asu ; (nanstd(nduration(zz))/sqrt(length(zz)))];
                end
             end
            %******** 
            iz = find( ~isnan(nuu));
            aa = [xuu(iz) ; flipud(xuu(iz))];
            bb = [(nuu(iz)+(2*nsu(iz))) ; flipud( nuu(iz)-(2*nsu(iz)))];
            fill(bb,aa,ncolo,'FaceAlpha',0.2,'Linestyle','none');
            plot(nuu(iz),xuu(iz),'k-','Color',ncolo,'Linewidth',2);
            %plot(nuu+(2*nsu),xuu,'k-','Color',ncolo,'Linewidth',1);
            %plot(nuu-(2*nsu),xuu,'k-','Color',ncolo,'Linewidth',1);
            %******
            iz = find( ~isnan(buu));
            aa = [xuu(iz) ; flipud(xuu(iz))];
            bb = [(buu(iz)+(2*bsu(iz))) ; flipud( buu(iz)-(2*bsu(iz)))];
            fill(bb,aa,bcolo,'FaceAlpha',0.2,'Linestyle','none');
            plot(buu,xuu,'k-','Color',bcolo,'Linewidth',2);
            % plot(buu+(2*bsu),xuu,'k-','Color',bcolo,'Linewidth',1);
            % plot(buu-(2*bsu),xuu,'k-','Color',bcolo,'Linewidth',1);
            %********
            %******
            if (0)
              iz = find( ~isnan(auu));
              aa = [xuu(iz) ; flipud(xuu(iz))];
              bb = [(auu(iz)+(2*asu(iz))) ; flipud( auu(iz)-(2*asu(iz)))];
              fill(bb,aa,[0,0,0],'FaceAlpha',0.2,'Linestyle','none');
              plot(auu,xuu,'k-','Linewidth',2);
            end
            %********
        end
        
        axis tight;
        V = axis;
        xmin = min(vx);
        xmax = max(vx);
        axis([xmin xmax V(3) V(4)]);
        V = axis;
        xmen = 0.5*(xmin+xmax);
        plot([xmen,xmen],[V(3),V(4)],'k-');
        plot([V(1),V(2)],[0,0],'r--');
        if (correlated)
          bz = find( ~isnan(nduration(narrow)) & ~isnan(ndepth(narrow)) );
         [r1,p1] = corr(nduration(narrow(bz)),ndepth(narrow(bz)),'type','Spearman');
          bz = find( ~isnan(nduration(broad)) & ~isnan(ndepth(broad)) );
          [r2,p2] = corr(nduration(broad(bz)),ndepth(broad(bz)),'type','Spearman');
          title(sprintf('Nw:%4.2f(%5.3f);Bd:%4.2f(%5.3f)',r1,p1,r2,p2));
        else
          title(sprintf('Narrow(%d);Broad(%d)',length(narrow),length(broad)));
        end
        subplot('Position',[xx 0.1 0.20 0.15]);
        nyx = hist(nduration(narrow),vx);
        nyx = nyx / sum(nyx);
        byx = hist(nduration(broad),vx);
        byx = byx / sum(byx);
        plot(vx,nyx,'k.-','Color',ncolo,'Linewidth',2); hold on;
        plot(vx,byx,'k.-','Color',bcolo,'Linewidth',2);
        axis tight;
        V = axis;
        axis([xmin xmax V(3) V(4)]);
        V = axis;
        plot([xmen,xmen],[V(3),V(4)],'k-');
        %***
        med1 = nanmedian(nduration(narrow));
        med2 = nanmedian(nduration(broad));
        mod1 = 100*(((1+med1)/(1-med1))-1);
        mod2 = 100*(((1+med2)/(1-med2))-1);
        p = ranksum(nduration(narrow),nduration(broad));
        if (correlated)
            title(sprintf('Nw:%4.1f,Bd:%4.1f(p=%5.4f)',mod1,mod2,p));
        end
return;

% non-parametric estimates of tuning parameters
function [ampa,ampu,basea,baseu,awid,uwid] = compute_NP_params(ratt,rutt)
          %**** assumes rank order around mean preference already done
%           basea = nanmean(ratt(1:4));
%           baseu = nanmean(rutt(1:4));
%           ampa = nanmean(ratt(13:16))-basea;
%           ampu = nanmean(rutt(13:16))-baseu;  
%           zaa = max((ratt-basea),0);
%           zaa = min((zaa/ampa),1);  % bound from 0 to 1
%           zuu = max((rutt-baseu),0);
%           zuu = min((zuu/ampu),1);  % bound from 0 to 1
%           awid = sum(zaa)/length(zaa);
%           uwid = sum(zuu)/length(zuu);
        
          basea = nanmean(ratt(1:1));
          baseu = nanmean(rutt(1:1));
          ampa = nanmean(ratt(16:16))-basea;
          ampu = nanmean(rutt(16:16))-baseu;  
          zaa = max((ratt-basea),0);
          zaa = min((zaa/ampa),1);  % bound from 0 to 1
          zuu = max((rutt-baseu),0);
          zuu = min((zuu/ampu),1);  % bound from 0 to 1
          %******* find mid-point in mass
          tot = sum(zaa);  % total
          otso = 0;
          tso = 0;
          mid = NaN;
          for i = 1:length(zaa)
              tso = (sum(zaa(1:i))/tot);
              if (tso > 0.5)
                 mid = ( i*(0.5-otso) + (i-1)*(tso-0.5) )/(tso-otso);
                 break;
              else
                 otso = tso; 
              end
          end
          awid = length(zaa)-mid;
          %********
          tot = sum(zuu);  % total
          otso = 0;
          tso = 0;
          mid = NaN;
          for i = 1:length(zuu)
              tso = sum(zuu(1:i))/tot;
              if (tso > 0.5)
                 mid = ( i*(0.5-otso) + (i-1)*(tso-0.5) )/(tso-otso);
                 break;
              else
                 otso = tso;
              end
          end
          uwid = length(zuu)-mid;   % from peak down to 0.5
          %*******  
return;


%***** scatter plot of AUC columns and stats
function plot_von_scatter(AUC,LockInt,WithLog,WithStd,Range,color)

   %going to need three scatters, one for each parameter
   
    ACounts = zeros(1,4);
    cmap = [[0.4,0,0.4];[0.7,0.3,0.7];[0,0.4,0];[0.3,0.7,0.3];[0.99,0.73,0.01];[0.96,0.815,0.42]];
    for k = 1:size(AUC,1)
        colo = cmap(color,:);
        %****** color code the dots      % for now, plot as one color
        if (AUC(k,1) == 1)
            if (AUC(k,2) == 1)
               % colo = cmap(1,:);
                ACounts(1) = ACounts(1)+1;
            else
               % colo = cmap(2,:);
                ACounts(2) = ACounts(2)+1;
            end
        else
            if (AUC(k,1) == 3)
               if (AUC(k,2) == 1)
                %   colo = cmap(3,:);
                   ACounts(3) = ACounts(3)+1;
               else
                 %  colo = cmap(4,:);
                   ACounts(4) = ACounts(4)+1;
               end
            else
               if (AUC(k,2) == 1)
                  % colo = cmap(5,:);
                   ACounts(1) = ACounts(1)+1;
               else
                  % colo = cmap(6,:);
                   ACounts(2) = ACounts(2)+1;
               end
            end
        end
        %*******
        if (WithLog)
           h = loglog(AUC(k,5),AUC(k,3),'k.'); hold on;
           set(h,'Markersize',12);
           set(h,'Color',colo);
           if (WithStd)
             h2 = loglog([(AUC(k,5)-AUC(k,6)),(AUC(k,5)+AUC(k,6))],[AUC(k,3),AUC(k,3)],'k-');
             set(h2,'Color',colo);
             h3 = loglog([AUC(k,5),AUC(k,5)],[(AUC(k,3)-AUC(k,4)),(AUC(k,3)+AUC(k,4))],'k-');
             set(h3,'Color',colo);
           end
        else
           h = plot(AUC(k,5),AUC(k,3),'k.'); hold on;
           set(h,'Markersize',12);
           set(h,'Color',colo);
           if (WithStd)
             h2 = plot([(AUC(k,5)-AUC(k,6)),(AUC(k,5)+AUC(k,6))],[AUC(k,3),AUC(k,3)],'k-');
             set(h2,'Color',colo);
             h3 = plot([AUC(k,5),AUC(k,5)],[(AUC(k,3)-AUC(k,4)),(AUC(k,3)+AUC(k,4))],'k-');
             set(h3,'Color',colo);
           end
        end
        %*********
    end    
    axis tight;
    V = axis;
    if isempty(Range)
      maxa = max(V(2),V(4));
      maxa = maxa*1.25;
      mina = min(V(1),V(3))*0.5;
    else
       maxa = Range(2);
       mina = Range(1);
    end
    axis([mina maxa mina maxa]);
    if (WithLog)
       loglog([mina,maxa],[mina,maxa],'k-'); % unity line
    else
       plot([mina,maxa],[mina,maxa],'k-'); % unity line
    end
    %****** comp basic stats
    ai = (AUC(:,3)-AUC(:,5))./(AUC(:,3)+AUC(:,5));
    medai = nanmedian(ai); % modo = (AUC(:,3) ./ AUC(:,5) );
    zz = find( ~isnan(ai) );  % problem is AI can be NaN if both zero rate
    medmodo = (1+medai)/(1-medai); % nanmedian(modo); 
    [p,h] = signrank(ai(zz)); % did not do pair obs otherwise % AUC(:,3),AUC(:,5));
    title(sprintf('Int(%4.2f,%4.2f) Mod:%5.3f (p=%6.4f)',...
                  LockInt(1),LockInt(2),medmodo,p));
return;

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

function [amed1,amed2] = compute_ai_print(ai1,ai2,topicname)

    amed1 = nanmedian(ai1);
    amed2 = nanmedian(ai2);
    med1 = (1+amed1)/(1-amed1);
    med2 = (1+amed2)/(1-amed2);
    med12 = mean([med1,med2]);
    [p1,h1] = signrank(ai1);
    [p2,h2] = signrank(ai2);
    [p12,h12] = ranksum(ai1,ai2);
    %******** disp stats in print
    disp([topicname,' comparison of AI indices']);
    disp(sprintf('MTC: med(%5.3f) ai(%5.3f) p(%6.4f)',med1,amed1,p1));
    disp(sprintf('MT: med(%5.3f) ai(%5.3f) p(%6.4f)',med2,amed2,p2));
    disp(sprintf('MTC vs M differ: med(%5.3f) animal differ p(%6.4f)',med12,p12));
    
return;
    
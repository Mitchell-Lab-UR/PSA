function includener = extra_inclusion_criteria(Info,StrictFit)
            
            %*** only include well fit von mises units
            POri = Info.OriTune.Afit.mu(4); % pref ori
            POriSD = Info.OriTune.Afit.sem(4); % sem
            %*** Then get att conditions
            APOri = Info.SacOn.OriTune{1}.AAfit.mu(4);
            APOriSD = Info.SacOn.OriTune{1}.AAfit.sem(4);
            IPOri = Info.SacOn.OriTune{2}.AAfit.mu(4);
            IPOriSD = Info.SacOn.OriTune{2}.AAfit.sem(4);
            %************
            if (abs(POri - APOri) > pi)
              APOri = APOri + sign(POri-APOri)*2*pi;
            end
            if (abs(POri - IPOri) > pi)
              IPOri = IPOri + sign(POri-IPOri)*2*pi;
            end
            %******
            POri = POri * (180/pi);
            APOri = APOri * (180/pi);
            IPOri = IPOri * (180/pi);
            POriSD = POriSD * (180/pi);
            APOriSD = APOriSD * (180/pi);
            IPOriSD = IPOriSD * (180/pi);
            % [POri,APOri,IPOri]
            % [POriSD,APOriSD,IPOriSD]
            %**********
            fname = Info.tagname;
            year = str2num(fname((end-4):end-3));
            month = str2num(fname((end-6):(end-5)));
            day = str2num(fname((end-8):(end-7)));
            %**** Only before July 2019 for Milo
%             if (0)
%                if ( (year > 19) || ((year == 19) && ((month >= 7)) ) ) 
%                   includener = 0;
%                   return;
%                end
%             end
%             if (1)
%                if (fname(1) == 'E')
%                  includener = 0;
%                  return;
%                end
%             end
            %*****
            if (StrictFit)  % selective for best tuned, many trials neurons
               ORI_THRESH = 360; 
               ORISD_THRESH = 360; 
               MIN_TRIALS = 150; 
            else  % include everything
               ORI_THRESH = 360; 
               ORISD_THRESH = 360;  % is this too harsh?
               MIN_TRIALS = 48;
            end
            %*** PFR inclusion requirement
            % if (Info.PFR_TargErr > 100)
            %    includener = 0;
            %    return;
            % end 
            
           %*** Kilosort inclusion
           if(1)
              if (Info.isolation(1) < 1)
                  includener = 0;
                  return;
              end
           end
           
           
            %*************
          %  includener = 1;
            if (1)
              ATrials = length(Info.AttList);
              ITrials = length(Info.IgnAList)+length(Info.IgnBList);
              if (ATrials < MIN_TRIALS) || (ITrials < MIN_TRIALS)
                  includener = 0;
                  return;
              end
            end
            %******* some high MU still left in dataset
            if (0)
              if (Info.BaseMu > 50) 
                 disp(sprintf('$$$$$$$$ Discarded unit %s for high MU',Info.tagname));
                 includener = 0;
                 return;               
              end
            end
            %******** more strict DSI criterion (insure pref ori strong)
            if (0)
              LowDSI = Info.OriTune.DSI - (2 * Info.OriTune.DSI_SEM);
              if (LowDSI < 0.1)
                 includener = 0;
                 return;
              end
            end
            %***** must match preferred direction, and must had good fit
            includener = 0;
            if ( abs(POri-APOri) < ORI_THRESH) && (abs(POri-IPOri) < ORI_THRESH)
                if (APOriSD < ORISD_THRESH) && (IPOriSD < ORISD_THRESH) 
                   includener = 1;
                end
            end
            %********
            if (0)  % show your plot to confirm it is good
              ovals = Info.SacOn.OriTune{1}.OriVals';
              araw = Info.SacOn.OriTune{1}.Mu;
              afit = Info.SacOn.OriTune{1}.AAtune.mu;
              iraw = Info.SacOn.OriTune{2}.Mu;
              ifit = Info.SacOn.OriTune{2}.AAtune.mu;
            
              figure(10); hold off;
              plot(ovals,araw,'ro'); hold on;
              plot(ovals,afit,'r-');
              plot(ovals,iraw,'bo'); hold on;
              plot(ovals,ifit,'b-');
              title(sprintf('Included %d',includener));
            
              input('check');
            end
            
 return;
            
          
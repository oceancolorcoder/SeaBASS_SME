% Input data from make_database_seabass.m
% Add AVW and QWIP if not already included in Derived Products
% Display spectral (ir)radiance and reflectance with ancillary data
% dashboard

% Screen on QWIP,RelAz,Wind,SZA, and visual inspection
%
% D. Aurin NASA/GSFC June 2024

%% Setup
wipe

% Set to true unless building .env.all for NOMAD/SeaBASS
ancillary.validation = 1;

ancillary.cruise = 'UMCES_Missouri_Reservoirs';
ancillary.SBA = 1;

clobber = 1;    % Re-evaluate thresholds and save over old
plotQWIP = 0;   % Plot QWIP (x2)
plotTimelineSpectra = 1;  % Plot timeline and spectral plot of flagged spectra
manualSelection = 1;    % Manually select/flag additional spectra

thresholds.negRrs = [380 680]; % Spectral range of negatives to eliminate from all sets
thresholds.relAz = [87 138]; % M99, Z17, IOCCG
thresholds.sza = [18 62]; % e.g. 20: Zhang 2017, depends on wind, e.g. 60:Brewin 2016
thresholds.wind = 10; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,
thresholds.qwip = 0.2; % Dierssen et al. 2022
thresholds.qa = 0.2; % This is more experimental. Monitor it, but don't filter it.
thresholds.cloud = [20 80]; % Clear and fully overcast should be okay. 20% - 80% are likely to be the worst. This is experimental.
thresholds.cloudIndexes = [0.05 0.3]; % From Ruddick et al. 2006 based on M99 models, where <0.05 is clear, >0.3 is fully overcast
if ancillary.validation
    % Thresholds for ancillary.validation
    thresholds.relAz = [89 136]; % M99, Z17, IOCCG
    thresholds.wind = 7; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,
    thresholds.qwip = 0.17; % Experimental
end

%% Load and overview
if clobber
    load(sprintf('dat/%s.mat',ancillary.cruise)) % dBase from make_awr_seabass.m
    AWR.dateTime = [dBase.datetime];
    AWR.nSpectra = length(dBase);

    try
        % Potential of wave mismatch depending on glint processing
        AWR.wave = dBase(1).wavelength;
        AWR.Rrs = vertcat(dBase.rrs);
        AWR.Es = vertcat(dBase.es);
        AWR.Lw = vertcat(dBase.lw);
        if ~ancillary.SBA
            % Have not encountered this yet
            AWR.Li = vertcat(dBase.li);            
        end
        if sum(contains(fieldnames(dBase),'rrs_sd')) > 0
            Rrs_sd = vertcat(dBase.rrs_sd);
            % Reduce resolution for clarity
            AWR.wave_sd = min(AWR.wave):10:max(AWR.wave);
            AWR.Rrs_sub = interp1(AWR.wave,AWR.Rrs',AWR.wave_sd)';
            AWR.Rrs_sd = interp1(AWR.wave,AWR.Rrs_sd',AWR.wave_sd)';
        end
        if sum(contains(fieldnames(dBase),'rrs_unc')) > 0
            % Note we are re-using _sd in the variable for uncertainty
            % (just for plots)
            Rrs_sd = vertcat(dBase.rrs_unc);
            AWR.wave_sd = min(AWR.wave):10:max(AWR.wave);
            AWR.Rrs_sub = interp1(AWR.wave,AWR.Rrs',AWR.wave_sd)';
            AWR.Rrs_sd = interp1(AWR.wave,AWR.Rrs_sd',AWR.wave_sd)';
        end
        clear Rrs_sd

    catch
        disp('Wavelength cannot be concatonated. Interpolating to shorter wave range')
        % Cross this bridge when/if it arises
    end
    %% Apply the QWIP and AVW, regardless of whether they are provided
    [ancillary.qwip,QCI,ancillary.avw] = AVW_QWIP_2D_fun(AWR.Rrs,AWR.wave,'none','none');
    if plotQWIP
        minMaxPlot = [440 590];
        fh12 = QWIP_figure_fun(AWR.Rrs,AWR.wave,ancillary.avw,ancillary.cruise,minMaxPlot);
        exportgraphics(fh12(1),sprintf('plt/%s_QWIP.png',ancillary.cruise))
        exportgraphics(fh12(2),sprintf('plt/%s_QWIP_Hist.png',ancillary.cruise))
    end

    %% Populate threshold variables
    if sum(contains(fieldnames(dBase),'chlor_a')) > 0
        ancillary.chl = [dBase.chlor_a];
    else
        % Calculate chl
        ancillary.chl=nan(1,AWR.nSpectra);
    end
    if sum(contains(fieldnames(dBase),'QA_score')) > 0
        ancillary.qa = [dBase.QA_score];
    else
        % Calculate QA_score
        ancillary.qa=nan(1,AWR.nSpectra);
    end
    if sum(contains(fieldnames(dBase),'cloud')) > 0
        ancillary.cloud = [dBase.cloud]; % Percent
    else
        ancillary.cloud = nan(1,AWR.nSpectra);
    end
    if ~ancillary.SBA
        % Fill in with cloud index where empty
        whrNaN = isnan(ancillary.cloud)';
        [wv,iwv] = find_nearest(750,AWR.wave);
        li750 = AWR.Li(:,iwv);
        es750 = AWR.Es(:,iwv);
        ancillary.cloud(whrNaN & (li750./es750 < thresholds.cloudIndexes(1))) = 0; % Clear
        ancillary.cloud(whrNaN & ( li750./es750 >= thresholds.cloudIndexes(1)) & li750./es750 < thresholds.cloudIndexes(2) ) = 50; % Partly cloudy
        ancillary.cloud(whrNaN & (li750./es750 > thresholds.cloudIndexes(2))) = 100; % Fully overcast (no good for validation, but okay for AWR)
    end

    if sum(contains(fieldnames(dBase),'wind')) > 0
        ancillary.wind = [dBase.wind];
    end
    if sum(contains(fieldnames(dBase),'SZA')) > 0
        ancillary.sza = [dBase.SZA];
    else
        disp('Calculating SZA')        
        % sun_position is not vector-ready
        for i=1:length(AWR.dateTime)
            [time.year,time.month,time.day,time.hour,time.min,time.sec] = datevec(AWR.dateTime(i));            
            time.UTC = 0*second(AWR.dateTime(i));
            location.latitude = [dBase(i).latitude];
            location.longitude = [dBase(i).longitude];

            sun = sun_position(time, location);
            ancillary.sza(i) = sun.zenith;
        end
    end
    if ~ancillary.SBA
        if sum(contains(fieldnames(SASData.Ancillary),'REL_AZ')) > 0
            ancillary.relAz = [dBase.REL_AZ];
        end
    end

    flags = writeFlags(ancillary,thresholds, AWR);
else
    load(sprintf('dat/%s_flags.mat',ancillary.cruise))
    AWR.nSpectra = length(AWR.Rrs); % check this
end

%% Show spectra with filtered data
if plotTimelineSpectra
    handles = plotFlags(AWR,ancillary,flags);


    %% Manual screening of spectra
    if manualSelection
        if ancillary.validation
            fprintf('Screen for validation\n')
        else
            fprintf('Screen for NOMAD\n')
        end

        input('Continue now? (enter)');
        % close all
        fh5 = figure;
        set(fh5,'position',[1926         381        1196         979])
        ah1 = axes;
        plot(AWR.wave,AWR.Rrs,'k')
        hold on
        flagSpectra(ah1,AWR.wave,AWR.Rrs,flags,0)
        ylabel('R_{rs} [sr^{-1}]')

        disp('Manual Screening')
        disp('Zoom to spectrum of interest and hit Continue')
        disp('Left mouse to flag spectrum. Continue to save and move on. Try again to ignore last.')
        disp('Middle mouse to exit and save.')
        flag = zeros(1,AWR.nSpectra);
        for tries=1:10
            h1 = uicontrol(fh5,'Style', 'pushbutton', 'String', 'Continue',...
                'Position', [5 100 50 25], 'Callback', 'uiresume');
            uicontrol(h1)
            uiwait(fh5)

            [x,y] = ginput(1);

            plot(x,y,'k*')
            h2 = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
                'Position', [5 100 50 25], 'Callback', 'butt=1; uiresume');
            h3 = uicontrol('Style', 'pushbutton', 'String', 'Try Again',...
                'Position', [5 50 50 25], 'Callback','butt=0; uiresume');
            h4 = uicontrol('Style', 'pushbutton', 'String', 'Exit',...
                'Position', [5 5 50 25], 'Callback','butt=3; uiresume');
            uicontrol(h2)
            uicontrol(h3)
            uicontrol(h4)
            uiwait
            if butt == 0
                continue
            elseif butt == 1
                [wv,windex] = find_nearest(x,AWR.wave);
                RrsX = AWR.Rrs(:,windex);
                [rrs,Rindex] = find_nearest(y,RrsX);
                plot(AWR.wave,AWR.Rrs(Rindex,:),'k','LineWidth',3)
                flag(Rindex) = 1;
                fprintf('Index of selected spectrum: %d\n',Rindex)
                disp(flag(Rindex))
                % disp([x,y])
                % continue
            elseif butt == 3
                break
            end
        end

        if size(flag,1)~=1
            flags.Manual = logical(flag');
        else
            flags.Manual = logical(flag);
        end
        %% Apply flags
        if ~ancillary.SBA
            flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.RelAz] | [flags.QWIP] |...
                [flags.Manual] | [flags.negRrs];
            set(handles.th8,'String',sprintf('Remaining: %d of %d',sum(~flag),AWR.nSpectra));
        else
            flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] |  [flags.QWIP] |...
            [flags.Manual] | [flags.negRrs];
            set(handles.th8,'String',sprintf('Remaining: %d of %d',sum(~flag),AWR.nSpectra));
        end

        flagSpectra(handles.ax11,AWR.wave,AWR.Rrs,flags,1)
        flagSpectra(handles.ax12,AWR.wave,AWR.Es,flags,0)
        if ~ancillary.SBA
            flagSpectra(ax13,AWR.wave,AWR.Li,flags,0)
        end
        flagSpectra(handles.ax14,AWR.wave,AWR.Lw,flags,0)

        if ancillary.validation
            exportgraphics(handles.fh3,sprintf('plt/%s_spec.png',ancillary.cruise))
        else
            exportgraphics(handles.fh3,sprintf('plt/%s_all_spec.png',ancillary.cruise))
        end
                
        %% Write CSV file for awr2env.py
        % In effect, we will have 2 files. One will have 0 or 1, the other 0 or 2.
        % Round to the nearest second
        dateTime = dateshift(AWR.dateTime,'start','minute') + seconds(round(second(AWR.dateTime)));
        [Yr,Mon,Day,Hr,Min,Sec] = datevec(dateTime');

        if ancillary.validation
            csvOutFile = sprintf('dat/%s_flags.csv',ancillary.cruise);
            FLAG = 2*int8(~flag'); % 0, 1, or 2 for reject, seabass-only, validation
        else
            csvOutFile = sprintf('dat/%s_all_flags.csv',ancillary.cruise);
            FLAG = int8(~flag'); % 0, 1, or 2 for reject, seabass-only, validation
        end
        T = table(Yr,Mon,Day,Hr,Min,Sec,FLAG);
        writetable(T,csvOutFile)
    end
end

% Resave with manual selections
if clobber
    save(sprintf('dat/%s_flags.mat',ancillary.cruise),"AWR","ancillary","flags")
end


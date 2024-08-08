% Input data from make_database_seabass.m
% Add AVW and QWIP if not already included in Derived Products
% Display spectral (ir)radiance and reflectance with ancillary data
% dashboard

% Screen on QWIP,RelAz,Wind,SZA, and visual inspection
%
% D. Aurin NASA/GSFC June 2024
%
%   Thresholds provided by the set_thresholds function

wipe
[fonts,projPath] = machine_prefs;
%% Manual Setup

% Set to true unless building .env.all for NOMAD/SeaBASS
ancillary.validation = 0;

% ancillary.cruise = 'BIOSCAPE_COASTAL_CARBON_Walker_Bay'; % Kyle Turner/Maria Tzortziou, BIOSCAPE (S. Africa)
% ancillary.cruise = 'BIOSCAPE_COASTAL_CARBON_St_Helena_Bay_2023'; % Kyle Turner/Maria Tzortziou, BIOSCAPE (S. Africa)
% ancillary.cruise = 'CCNY_tzortziou_ARCTICCC_Norton_Sound_2022_AWRrrs';
% ancillary.cruise = 'CCNY_tzortziou_ARCTICCC_Norton_Sound_2023_AWRrrs';
% ancillary.cruise = 'ArcticCC_Norton_Sound_2022';
% ancillary.cruise = 'ArcticCC_Alakanuk_2022';
% ancillary.cruise = 'ArcticCC_Alakanuk_2023';
% ancillary.cruise = 'ArcticCC_Norton_Sound_2023';
ancillary.cruise = 'VIIRS_VALIDATION_viirs_2021_gunter';

ancillary.SBA = 0;

% relAz = (135+90)/2;   % Not provided per cast, but in this range
relAz = 90;             % Not provided per cast, but reported as about +/-90
% ancillary.cruise = 'UMCES_Missouri_Reservoirs';
% ancillary.SBA = 1;

clobber = 1;                % Re-evaluate thresholds and save over old
plotQWIP = 0;               % Plot QWIP (x2)
plotTimelineSpectra = 1;    % Req'd for manualSelection. Plot timeline and spectral plot of flagged spectra
manualSelection = 1;        % Manually select/flag additional spectra

%% Auto Setup
thresholds = set_thresholds(ancillary.validation);
SMEPath = fullfile(projPath,'SeaBASS','JIRA_Tickets',ancillary.cruise);
plotPath = fullfile(SMEPath,'Plots');
if ~isfolder(plotPath)
    mkdir(plotPath)
end

%% Load and overview
if clobber
    load(sprintf('dat/%s.mat',ancillary.cruise)) % dBase from make_awr_seabass.m
    AWR.dateTime = [dBase.datetime];
    AWR.nSpectra = length(dBase);

    try
        % Potential of wave mismatch depending on glint processing
        AWR.wave = dBase(1).wavelength;
    catch
        disp('***** Wavelength cannot be concatonated. Interpolating to shorter wave range ****')
        % Cross this bridge when/if it arises
    end
    AWR.Rrs = vertcat(dBase.rrs);
    AWR.Es = vertcat(dBase.es);
    if isfield(dBase,'lw')
        AWR.Lw = vertcat(dBase.lw);
    end
    if isfield(dBase,'lt')
        AWR.Lt = vertcat(dBase.lt);
    end
    if isfield(dBase,'li')
        AWR.Li = vertcat(dBase.li);
    end
    AWR.wave_sd = min(AWR.wave):10:max(AWR.wave);
    if isfield(dBase,'rrs_sd')
        Rrs_sd = vertcat(dBase.rrs_sd);
        % Reduce resolution for clarity        
        AWR.Rrs_sub = interp1(AWR.wave,AWR.Rrs',AWR.wave_sd)';
        AWR.Rrs_sd = interp1(AWR.wave,Rrs_sd',AWR.wave_sd)';
    elseif isfield(dBase,'rrs_unc')
        % Note we are re-using _sd in the variable for uncertainty
        % (just for plots)
        Rrs_sd = vertcat(dBase.rrs_unc);
        AWR.Rrs_sub = interp1(AWR.wave,AWR.Rrs',AWR.wave_sd)';
        AWR.Rrs_sd = interp1(AWR.wave,Rrs_sd',AWR.wave_sd)';
    end
    if isfield(dBase,'es_sd')
        Es_sd = vertcat(dBase.es_sd);
        % Reduce resolution for clarity
        AWR.Es_sub = interp1(AWR.wave,AWR.Es',AWR.wave_sd)';
        AWR.Es_sd = interp1(AWR.wave,Es_sd',AWR.wave_sd)';
    elseif isfield(dBase,'es_unc')
        % Note we are re-using _sd in the variable for uncertainty
        % (just for plots)
        Es_sd = vertcat(dBase.es_unc);
        AWR.Es_sub = interp1(AWR.wave,AWR.Es',AWR.wave_sd)';
        AWR.Es_sd = interp1(AWR.wave,Es_sd',AWR.wave_sd)';
    end
    if isfield(dBase,'lt_sd')
        Lt_sd = vertcat(dBase.lt_sd);
        % Reduce resolution for clarity
        AWR.Lt_sub = interp1(AWR.wave,AWR.Lt',AWR.wave_sd)';
        AWR.Lt_sd = interp1(AWR.wave,Lt_sd',AWR.wave_sd)';
    elseif isfield(dBase,'lt_unc')
        % Note we are re-using _sd in the variable for uncertainty
        % (just for plots)
        Lt_sd = vertcat(dBase.lt_unc);
        AWR.Lt_sub = interp1(AWR.wave,AWR.Lt',AWR.wave_sd)';
        AWR.Lt_sd = interp1(AWR.wave,Lt_sd',AWR.wave_sd)';
    end
    if isfield(dBase,'lw_sd')
        Lw_sd = vertcat(dBase.lw_sd);
        % Reduce resolution for clarity
        AWR.Lw_sub = interp1(AWR.wave,AWR.Lw',AWR.wave_sd)';
        AWR.Lw_sd = interp1(AWR.wave,Lw_sd',AWR.wave_sd)';
    elseif isfield(dBase,'lw_unc')
        % Note we are re-using _sd in the variable for uncertainty
        % (just for plots)
        Lw_sd = vertcat(dBase.lw_unc);
        AWR.Lw_sub = interp1(AWR.wave,AWR.Lw',AWR.wave_sd)';
        AWR.Lw_sd = interp1(AWR.wave,Lw_sd',AWR.wave_sd)';
    end
    if isfield(dBase,'li_sd')
        Li_sd = vertcat(dBase.li_sd);
        % Reduce resolution for clarity
        AWR.Li_sub = interp1(AWR.wave,AWR.Li',AWR.wave_sd)';
        AWR.Li_sd = interp1(AWR.wave,Li_sd',AWR.wave_sd)';
    elseif isfield(dBase,'li_unc')
        % Note we are re-using _sd in the variable for uncertainty
        % (just for plots)
        Li_sd = vertcat(dBase.li_unc);
        AWR.Li_sub = interp1(AWR.wave,AWR.Li',AWR.wave_sd)';
        AWR.Li_sd = interp1(AWR.wave,Li_sd',AWR.wave_sd)';
    end
    clear Rrs_sd

    
    %% Apply the QWIP and AVW, regardless of whether they are provided
    [ancillary.qwip,QCI,ancillary.avw] = AVW_QWIP_2D_fun(AWR.Rrs,AWR.wave,'none','none');
    if plotQWIP
        minMaxPlot = [440 590];
        fh12 = QWIP_figure_fun(AWR.Rrs,AWR.wave,ancillary.avw,ancillary.cruise,minMaxPlot);
        exportgraphics(fh12(1),sprintf('%s/%s_QWIP.png',plotPath,ancillary.cruise))
        exportgraphics(fh12(2),sprintf('%s/%s_QWIP_Hist.png',plotPath,ancillary.cruise))
    end

    %% Populate threshold variables
    if isfield(dBase,'chlor_a')
        ancillary.chl = [dBase.chlor_a];
    else
        % Calculate chl
        ancillary.chl=nan(1,AWR.nSpectra);
    end
    if isfield(dBase,'QA_score')
        ancillary.qa = [dBase.QA_score];
    else
        % Calculate QA_score
        ancillary.qa=nan(1,AWR.nSpectra);
    end
    if isfield(dBase,'cloud')
        ancillary.cloud = [dBase.cloud]; % Percent
    else
        ancillary.cloud = nan(1,AWR.nSpectra);
    end
    if ~ancillary.SBA
        if ~isfield(AWR,'Li')
            disp('****Check that this is not SBA****')
            return
        end
        % Fill in with cloud index where empty
        whrNaN = isnan(ancillary.cloud)';
        [wv,iwv] = find_nearest(750,AWR.wave);
        li750 = AWR.Li(:,iwv);
        es750 = AWR.Es(:,iwv);
        ancillary.cloud(whrNaN & (li750./es750 < thresholds.cloudIndexes(1))) = 0; % Clear
        ancillary.cloud(whrNaN & ( li750./es750 >= thresholds.cloudIndexes(1)) & li750./es750 < thresholds.cloudIndexes(2) ) = 50; % Partly cloudy
        ancillary.cloud(whrNaN & (li750./es750 > thresholds.cloudIndexes(2))) = 100; % Fully overcast (no good for validation, but okay for AWR)
    end

    if isfield(dBase,'wind')
        ancillary.wind = [dBase.wind];
    end
    if isfield(dBase,'SZA')
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
        if exist('relAz','var')
            ancillary.relAz = repmat(relAz,size(ancillary.sza));
        end
    end

    flags = writeFlags(ancillary,thresholds, AWR);
else
    load(sprintf('dat/%s_flags.mat',ancillary.cruise))
    AWR.nSpectra = length(AWR.Rrs); % check this
end

%% Show spectra with filtered data
if plotTimelineSpectra    
    handles = plotFlags(AWR,ancillary,flags,plotPath);

    handles.fh5 = figure('position',[1926         381        1196         979]);
    handles.ah1 = axes;    
    plot(AWR.wave,AWR.Rrs,'k')
    hold on
    if isfield(AWR,'Rrs_sd')
        errorbar(AWR.wave_sd,AWR.Rrs_sub,AWR.Rrs_sd,'color',[0.8 0.8 0.8])
    end
    grid on
    flagSpectra(handles.ah1,AWR.wave,AWR.Rrs,flags,0)
    ylabel('R_{rs} [sr^{-1}]')

    %% Manual screening of spectra
    if manualSelection
        flags.Manual = manualFlag(ancillary,handles,AWR,flags,plotPath);
        set(handles.th0, 'String', sprintf('Manual: %d',sum(flags.Manual)),'Units','normalized');

        if ancillary.validation
            exportgraphics(handles.fh3,sprintf('%s/%s_spec.png',plotPath,ancillary.cruise))
        else
            exportgraphics(handles.fh3,sprintf('%s/%s_all_spec.png',plotPath,ancillary.cruise))
        end
    end
end

% Resave with manual selections
if clobber
    save(sprintf('dat/%s_flags.mat',ancillary.cruise),"AWR","ancillary","flags")
end


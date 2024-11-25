% Input data from make_database_hypercp.m
% Add AVW and QWIP if not already included in Derived Products
% Display spectral (ir)radiance and reflectance with ancillary data
% dashboard
%
% Screen on QWIP,RelAz,Wind,SZA, and visual inspection
%
%   NOTE: See ./sub/set_thresholds.m for QC thresholds and citations
%
%   Inputs:
%       dat/[cruisename].mat from make_awr_hypercp.m (SASData structure)
%
%   Outputs:
%       dat/[cruisename]_flags.mat (select variables and QA/QC flags)
%       plt/[cruisename]_AllSpec.png (Rrs, Es, Li, Lw with flags shown)
%       plt/[cruisename]_QWIP.png
%       plt/[cruisename]_QWIP_Hist.png
%       plt/[cruisename]_timeline.png (time series of select variables used
%       in QAQC)
%       dat/[cruisename]_all_flags.mat
%       dat/[cruisename]_flags.mat
%       dat/[cruisename]_all_flags.csv
%       dat/[cruisename]_flags.csv
%
% D. Aurin NASA/GSFC November 2024

% path(path,'./sub') % <-- uncomment if you are not me
wipe
[fonts,projPath] = machine_prefs;
%% Manual Setup

% Set to true unless building .env.all for NOMAD/SeaBASS
ancillary.validation = 1;
ancillary.cruise = 'NORTHERN_INDIAN_OCEAN_EKAMSAT-EKAMSAT-2024-Bay-of-Bengal';
% ancillary.cruise = 'KORUS_KR_2016_RV_Onnuri_HyperSAS';
% ancillary.cruise = 'EXPORTS_EXPORTSNP_Mannino_AOP_HyperSAS_R0';
% ancillary.cruise = 'EXPORTSNA_NASA';
% ancillary.cruise = 'EXPORTSNA_Boss';
% ancillary.cruise = 'Brewin_Superyacht_Science_2019-2020';
% ancillary.cruise = 'Brewin_Superyacht_Science_2018';

clobber = 1;    % Re-evaluate thresholds and save over old
plotQWIP = 1;   % Plot QWIP (x2)
plotTimelineSpectra = 1;  % Plot timeline and spectral plot of flagged spectra
manualSelection = 1;    % Manually select/flag additional spectra

%% Auto Setup
thresholds = set_thresholds(ancillary.validation);
ancillary.SBA = 0; % Set to 1 for SBA data (used by subroutines for review_awr_seabass, not here)
SMEPath = fullfile(projPath,'SeaBASS','JIRA_Tickets',ancillary.cruise);
plotPath = fullfile(SMEPath,'Plots');
if ~isfolder(plotPath)
    mkdir(plotPath)
end

%% Load and overview
if clobber
    load(sprintf('dat/%s.mat',ancillary.cruise)) % SASData from make_database_hypercp.m
    AWR.dateTime = datetime([SASData.Ancillary.datenum],'ConvertFrom','datenum','TimeZone','UTC');
    AWR.nSpectra = size([SASData.Ancillary.datenum],2);

    try
        % Potential of wave mismatch depending on glint processing
        AWR.wave = vertcat(SASData.Reflectance(1).Rrs_HYPER_wavelength);
        AWR.rrs = vertcat(SASData.Reflectance.Rrs_HYPER);
        AWR.es = vertcat(SASData.Irradiance.ES_HYPER);
        AWR.li = vertcat(SASData.Radiance.LI_HYPER);
        AWR.lw = vertcat(SASData.Radiance.LW_HYPER);
    catch
        disp('Wavelength cannot be concatonated. Interpolating to shorter wave range')
        % In this case, interpolation of data to the shorter wavelength vector
        % will be required
        % There has to be a better way...
        for i=AWR.nSpectra:-1:1
            nBands(i) = length(SASData.Reflectance(i).Rrs_HYPER_wavelength);
        end
        iShort = find(nBands == min(nBands)); iShort = iShort(1);
        AWR.wave = SASData.Reflectance(iShort).Rrs_HYPER_wavelength;
        for i=AWR.nSpectra:-1:1
            AWR.rrs(i,:) = interp1(SASData.Reflectance(i).Rrs_HYPER_wavelength, ...
                SASData.Reflectance(i).Rrs_HYPER, AWR.wave);
            AWR.es(i,:) = interp1(SASData.Irradiance(i).ES_HYPER_wavelength, ...
                SASData.Irradiance(i).ES_HYPER, AWR.wave);
            AWR.li(i,:) = interp1(SASData.Radiance(i).LI_HYPER_wavelength, ...
                SASData.Radiance(i).LI_HYPER, AWR.wave);
            AWR.lw(i,:) = interp1(SASData.Radiance(i).LW_HYPER_wavelength, ...
                SASData.Radiance(i).LW_HYPER, AWR.wave);
        end
    end
    %% Apply the QWIP and AVW, regardless of whether they are provided
    [ancillary.qwip,QCI,ancillary.avw] = AVW_QWIP_2D_fun(AWR.rrs,AWR.wave,'none','none');
    if plotQWIP
        minMaxPlot = [440 590];
        fh12 = QWIP_figure_fun(AWR.rrs,AWR.wave,ancillary.avw,ancillary.cruise,minMaxPlot);
        exportgraphics(fh12(1),sprintf('%s/%s_QWIP.png',plotPath,ancillary.cruise))
        exportgraphics(fh12(2),sprintf('%s/%s_QWIP_Hist.png',plotPath,ancillary.cruise))
    end

    %% Populate threshold variables
    if sum(contains(fieldnames(SASData),'Derived_Products')) > 0
        if sum(contains(fieldnames(SASData.Derived_Products),'chlor_a')) > 0
            ancillary.chl = [SASData.Derived_Products.chlor_a];
        else
            % Calculate chl
            ancillary.chl=nan;
        end
    else
        % Calculate chl
        ancillary.chl=nan;
    end
    if sum(contains(fieldnames(SASData),'Derived_Products')) > 0
        if sum(contains(fieldnames(SASData.Derived_Products),'QA_score')) > 0
            ancillary.qa = [SASData.Derived_Products.QA_score];
        else
            % Calculate QA_score
            ancillary.qa=nan;
        end
    else
        % Calculate QA_score
        ancillary.qa=nan;
    end
    if sum(contains(fieldnames(SASData.Ancillary),'CLOUD')) > 0
        ancillary.cloud = [SASData.Ancillary.CLOUD]; % Percent

        % Fix bad interpolation in EXPORTSNA Ancillary file
        ancillary.cloud(ancillary.cloud<0) = missing;
    else
        ancillary.cloud = nan(1,AWR.nSpectra);
    end
    % Fill in with cloud index where empty
    whrNaN = isnan(ancillary.cloud)';
    [wv,iwv] = find_nearest(750,AWR.wave);
    li750 = AWR.li(:,iwv);
    es750 = AWR.es(:,iwv);
    ancillary.cloud(whrNaN & (li750./es750 < thresholds.cloudIndexes(1))) = 0; % Clear
    ancillary.cloud(whrNaN & ( li750./es750 >= thresholds.cloudIndexes(1)) & li750./es750 < thresholds.cloudIndexes(2) ) = 50; % Partly cloudy
    ancillary.cloud(whrNaN & (li750./es750 > thresholds.cloudIndexes(2))) = 100; % Fully overcast (no good for validation, but okay for AWR)

    if sum(contains(fieldnames(SASData.Ancillary),'WINDSPEED')) > 0
        ancillary.wind = [SASData.Ancillary.WINDSPEED];
    end
    if sum(contains(fieldnames(SASData.Ancillary),'SZA')) > 0
        ancillary.sza = [SASData.Ancillary.SZA];
    end
    if sum(contains(fieldnames(SASData.Ancillary),'REL_AZ')) > 0
        ancillary.relAz = [SASData.Ancillary.REL_AZ];
    end

    flags = writeFlags(ancillary,thresholds, AWR);
else
    load(sprintf('dat/%s_flags.mat',ancillary.cruise))
    AWR.nSpectra = size(AWR.rrs,1);
end

%% Show spectra with filtered data
if plotTimelineSpectra
    handles = plotFlags(AWR,ancillary,flags,plotPath);

    handles.fh5 = figure('position',[1926         381        1196         979]);
    handles.ah1 = axes;
    plot(AWR.wave,AWR.rrs,'k')
    hold on
    grid on
    flagSpectra(handles.ah1,AWR.wave,AWR.rrs,flags,0)
    ylabel('R_{rs} [sr^{-1}]')

    %% Manual screening of spectra
    if manualSelection
        manualFlag(ancillary,handles,AWR,flags,plotPath);

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

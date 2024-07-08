% Input data from make_database_hypercp.m
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

% ancillary.cruise = 'EXPORTSNA_NASA';
ancillary.cruise = 'EXPORTSNA_Boss';
% ancillary.cruise = 'Brewin_Superyacht_Science_2019-2020';

ancillary.SBA = 0;

clobber = 1;    % Re-evaluate thresholds and save over old
plotQWIP = 1;   % Plot QWIP (x2)
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
    load(sprintf('dat/%s.mat',ancillary.cruise)) % SASData from make_database_hypercp.m
    AWR.dateTime = datetime([SASData.Ancillary.datenum],'ConvertFrom','datenum','TimeZone','UTC');
    AWR.nSpectra = size([SASData.Ancillary.datenum],2);

    try
        % Potential of wave mismatch depending on glint processing
        AWR.wave = vertcat(SASData.Reflectance(1).Rrs_HYPER_wavelength);
        AWR.Rrs = vertcat(SASData.Reflectance.Rrs_HYPER);
        AWR.Es = vertcat(SASData.Irradiance.ES_HYPER);
        AWR.Li = vertcat(SASData.Radiance.LI_HYPER);
        AWR.Lw = vertcat(SASData.Radiance.LW_HYPER);
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
            AWR.Rrs(i,:) = interp1(SASData.Reflectance(i).Rrs_HYPER_wavelength, ...
                SASData.Reflectance(i).Rrs_HYPER, AWR.wave);
            AWR.Es(i,:) = interp1(SASData.Irradiance(i).ES_HYPER_wavelength, ...
                SASData.Irradiance(i).ES_HYPER, AWR.wave);
            AWR.Li(i,:) = interp1(SASData.Radiance(i).LI_HYPER_wavelength, ...
                SASData.Radiance(i).LI_HYPER, AWR.wave);
            AWR.Lw(i,:) = interp1(SASData.Radiance(i).LW_HYPER_wavelength, ...
                SASData.Radiance(i).LW_HYPER, AWR.wave);
        end
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
    li750 = AWR.Li(:,iwv);
    es750 = AWR.Es(:,iwv);
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
    AWR.nSpectra = size(AWR.Rrs,1);
end

%% Show spectra with filtered data
if plotTimelineSpectra
    handles = plotFlags(AWR,ancillary,flags);    

    handles.fh5 = figure('position',[1926         381        1196         979]);
    handles.ah1 = axes;
    plot(AWR.wave,AWR.Rrs,'k')
    hold on
    grid on
    flagSpectra(handles.ah1,AWR.wave,AWR.Rrs,flags,0)
    ylabel('R_{rs} [sr^{-1}]')

    %% Manual screening of spectra
    if manualSelection
        manualFlag(ancillary,handles,AWR,flags)
    end
end

% Resave with manual selections
if clobber
    save(sprintf('dat/%s_flags.mat',ancillary.cruise),"AWR","ancillary","flags")
end

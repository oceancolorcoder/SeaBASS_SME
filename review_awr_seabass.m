% Input data from make_database_seabass.m
% Add AVW and QWIP if not already included in Derived Products
% Display spectral (ir)radiance and reflectance with ancillary data
% dashboard

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
%
%   Thresholds provided by the set_thresholds function

% path(path,'./sub') % <-- uncomment if you are not me
wipe
[fonts,projPath] = machine_prefs;   % <-- Set this
%% Manual Setup

% Set to true (1) unless building .env.all for NOMAD/SeaBASS (0)
ancillary.validation = 1;

% ancillary.cruise = 'BIOSCAPE_COASTAL_CARBON_Walker_Bay';          % Kyle Turner/Maria Tzortziou, BIOSCAPE (S. Africa)
% ancillary.cruise = 'BIOSCAPE_COASTAL_CARBON_St_Helena_Bay_2023';    % Kyle Turner/Maria Tzortziou, BIOSCAPE (S. Africa)
% ancillary.cruise = 'CCNY_tzortziou_ARCTICCC_Norton_Sound_2022_AWRrrs';
% ancillary.cruise = 'CCNY_tzortziou_ARCTICCC_Norton_Sound_2023_AWRrrs';
% ancillary.cruise = 'ArcticCC_Norton_Sound_2022';
% ancillary.cruise = 'ArcticCC_Alakanuk_2022';
% ancillary.cruise = 'ArcticCC_Alakanuk_2023';
% ancillary.cruise = 'ArcticCC_Norton_Sound_2023';
% ancillary.cruise = 'VIIRS_VALIDATION_viirs_2021_gunter';
% ancillary.cruise = 'VIIRS_VALIDATION_viirs_2022_sette';
% ancillary.cruise = 'VIIRS_VALIDATION_viirs_2023_shimada';
% ancillary.cruise = 'PVST_PRINGLS_PRINGLS_20240417';
% ancillary.cruise = 'PVST_PRINGLS_PRINGLS_20240513';
% ancillary.cruise = 'PVST_PRINGLS_PRINGLS_20240612';
% ancillary.cruise = 'UMCES_Missouri_Reservoirs';
% ancillary.cruise = 'NF2405_VIIRS';
% ancillary.cruise = 'CHESAPEAKE_BAY_HELICOPTER_Chesapeake_Bay_2022';
% ancillary.cruise = 'PVST_PRINGLS_PRINGLS_20240717';
ancillary.cruise = 'PVST_PRINGLS_PRINGLS_20240813';
% ancillary.cruise = 'PVST_PRINGLS_PRINGLS_20241003';
% ancillary.cruise = 'PVST_PRINGLS_PRINGLS_20240911';

SMEPath = fullfile(projPath,'SeaBASS','JIRA_Tickets',ancillary.cruise); % <-- Set this; used to write plots

ancillary.SBA = 0;                                          % <-- Set this (1 for SBA)
ancillary.skipLi = 0;                 % Used for cloud index. May not be available, even for non-SBA

% relAz = 90;               % If not provided per cast, but reported as +/-90   % <-- Check this
% relAz = 135;              % If not provided per cast, but reported as +/-135
% relAz = (135+90)/2;       % If not provided per cast, but in this range; Bad choice...     

clobber = 1;                % Re-evaluate thresholds and save over old
plotQWIP = 1;               % Plot QWIP (x2)
plotTimelineSpectra = 1;    % Req'd for manualSelection. Plot timeline and spectral plot of flagged spectra
manualSelection = 1;        % Manually select/flag additional spectra

%% Auto Setup
thresholds = set_thresholds(ancillary.validation);
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
        AWR.wave = dBase(1).wavelength;
    catch
        AWR.wave = dBase(1).rrs_wavelength;
    end
    
    AWR.wave_sd = min(AWR.wave):20:max(AWR.wave); % Subset for errorbars
    try
        % Potential of wave mismatch depending on glint processing
        AWR = vertcatFun(AWR,dBase);        
    catch
        disp('***** Wavelength cannot be concatonated. Interpolating to shorter wave range ****')
        AWR = interpFun(AWR,dBase);
    end    


    %% Apply the QWIP and AVW, regardless of whether they are provided
    [ancillary.qwip,QCI,ancillary.avw] = AVW_QWIP_2D_fun(AWR.rrs,AWR.wave,'none','none');
    if plotQWIP
        minMaxPlot = [440 590];
        
        % Figures 1 & 2
        fh12 = QWIP_figure_fun(AWR.rrs,AWR.wave,ancillary.avw,ancillary.cruise,minMaxPlot);
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
    if ~ancillary.SBA && ~ancillary.skipLi
        if ~isfield(AWR,'li')
            disp('****Check that this is not SBA****')
            return
        end
        % Fill in with cloud index where empty
        whrNaN = isnan(ancillary.cloud)';
        [wv,iwv] = find_nearest(750,AWR.wave);
        li750 = AWR.li(:,iwv);
        es750 = AWR.es(:,iwv);
        ancillary.cloud(whrNaN & (li750./es750 < thresholds.cloudIndexes(1))) = 0; % Clear
        ancillary.cloud(whrNaN & ( li750./es750 >= thresholds.cloudIndexes(1)) & li750./es750 < thresholds.cloudIndexes(2) ) = 50; % Partly cloudy
        ancillary.cloud(whrNaN & (li750./es750 > thresholds.cloudIndexes(2))) = 100; % Fully overcast (no good for validation, but okay for AWR)
    end

    if isfield(dBase,'wind')
        ancillary.wind = [dBase.wind];
    end
    if isfield(dBase,'sza')
        ancillary.sza = [dBase.sza];
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
    if isfield(dBase,'relaz')
        ancillary.relAz = [dBase.relaz];
    elseif isfield(dBase,'relAz')
        ancillary.relAz = [dBase.relAz];
    else
        if ~ancillary.SBA
            if exist('relAz','var')
                ancillary.relAz = repmat(relAz,size(ancillary.sza));
            end
        end
    end

    flags = writeFlags(ancillary,thresholds, AWR);
else
    load(sprintf('dat/%s_flags.mat',ancillary.cruise))
    AWR.nSpectra = length(AWR.rrs); % check this
end

%% Show spectra and timeline with filtered data
if plotTimelineSpectra
    
    % Figures 3 (spectral) & 4 (timeline)
    handles = plotFlags(AWR,ancillary,flags,plotPath);

    % Figure 5 (manual)
    handles.fh5 = figure('position',[1926         381        1196         979]);
    handles.ah1 = axes;
    plot(AWR.wave,AWR.rrs,'k')
    hold on
    if isfield(AWR,'rrs_sd')
        errorbar(AWR.wave_sd,AWR.rrs_sub,AWR.rrs_sd,'color',[0.8 0.8 0.8])
    end
    grid on
    flagSpectra(handles.ah1,AWR.wave,AWR.rrs,flags,0)
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


%%
function AWR = vertcatFun(AWR,dBase)

gudFields = ["rrs" "rrs_unc" "rrs_sd" "es" "es_unc" "es_sd" "lw" "lw_unc" "lw_sd" "lt" "lt_unc" "lt_sd" "li" "li_sd" "li_unc"];
dBaseFields = fieldnames(dBase);

for namei=1:length(dBaseFields)
    % Vectorized operations with vertcat
    fieldName = dBaseFields{namei};
    if any(strcmpi(fieldName, gudFields))
        % For uncertainty fields, interpolated to course spectral resolution
        if contains(fieldName,"_unc") || contains(fieldName,"_sd")
            AWR.(fieldName) = interp1(dBase.wavelength,vertcat(dBase.(fieldName))',AWR.wave_sd)';
        else
            % Otherwise full resolution and add low resolution for errorbar plots
            AWR.(fieldName) = vertcat(dBase.(fieldName));
            subField = [fieldName '_sub'];
            AWR.(subField) = interp1(AWR.wave,AWR.(subfield)',AWR.wave_sd)';            
        end
    end
end

end

%%
function AWR = interpFun(AWR,dBase)
% Automatically interpolates all to the first wavelength in dBase, if
% necessary

disp('Interpolating spectral results with differing wavebands')

wavelength = AWR.wave;
gudFields = ["rrs" "rrs_unc" "rrs_sd" "es" "es_unc" "es_sd" "lw" "lw_unc" "lw_sd" "lt" "lt_unc" "lt_sd" "li" "li_sd" "li_unc"];
dBaseFields = fieldnames(dBase);

for namei=1:length(dBaseFields)
    fieldName = dBaseFields{namei};
    if any(strcmpi(fieldName, gudFields))

        % Row-by-row interpolation
        for i=length(dBase):-1:1  
            try
                test = length(dBase(i).wavelength) ~= length(wavelength);
                waveOther = dBase(i).wavelength;
            catch
                fieldNameWv = sprintf('%s_wavelength',fieldName);
                fieldNameWv = strrep(fieldNameWv,'_sd_','_');
                fieldNameWv = strrep(fieldNameWv,'_unc_','_');
                test = length(dBase(i).(fieldNameWv)) ~= length(wavelength);
                waveOther = dBase(i).(fieldNameWv);
            end
            % For uncertainty fields, interpolated to coarse spectral resolution
            if contains(fieldName,"_unc") || contains(fieldName,"_sd")        
                AWR.(fieldName)(i,:) = interp1(waveOther, dBase(i).(fieldName), AWR.wave_sd);
            else
                % Otherwise full resolution
                
                if test
                    AWR.(fieldName)(i,:) = interp1(waveOther, dBase(i).(fieldName), wavelength);
                else
                    AWR.(fieldName)(i,:) = dBase(i).(fieldName);
                end
                % And add a low resolution version for errorbar plots
                subField = [fieldName '_sub'];
                AWR.(subField)(i,:) = interp1(waveOther, dBase(i).(fieldName), AWR.wave_sd);
            end

        end
    end
end



end
% Using readsb.m (in MATLAB_embedded, updated No 2023), pull in all the
% files obtained using run_get_files.bsh.
% 
% Inputs:
%        cruisename
%        SeaBASS files from .lis and run_get_files.bsh
%
% Output:
%   dBase structure with relevant fields (e.g., datetime, Rrs, etc.)
%
% D. Aurin NASA/GSFC November 2024

% path(path,'./sub')          % <-- uncomment if you are not me
wipe
%% Setup

% Cruise Name
%   The folder of SeaBASS files:

% cruise = 'viirs_2019_foster';
% cruise = 'RSWQ_2023'; % single spectrum per file
% cruise = 'JackBlanton'; % Rivero-Calle; returned to PI
% cruise = 'Belgium_2021'; % Twardowski; Not AWR. Returned to Mike.
% cruise = 'UMCES_Missouri_Reservoirs'; % SBA Lorena Silva/Greg Silsbe
% cruise = 'Brewin_Superyacht_Science_2018';
% cruise = 'BIOSCAPE_COASTAL_CARBON_Walker_Bay'; % Kyle Turner/Maria Tzortziou, BIOSCAPE (S. Africa)
% cruise = 'BIOSCAPE_COASTAL_CARBON_St_Helena_Bay_2023';
% cruise = 'CCNY_tzortziou_ARCTICCC_Norton_Sound_2022_AWRrrs';
% cruise = 'CCNY_tzortziou_ARCTICCC_Norton_Sound_2023_AWRrrs';
% cruise = 'ArcticCC_Norton_Sound_2022';
% cruise = 'ArcticCC_Alakanuk_2022';
% cruise = 'ArcticCC_Alakanuk_2023';
% cruise = 'ArcticCC_Norton_Sound_2023';
% cruise = 'VIIRS_VALIDATION_viirs_2021_gunter';
% cruise = 'VIIRS_VALIDATION_viirs_2023_shimada';
% cruise = 'VIIRS_VALIDATION_viirs_2022_sette';
% cruise = 'PVST_PRINGLS_PRINGLS_20240417';
% cruise = 'PVST_PRINGLS_PRINGLS_20240513';
% cruise = 'PVST_PRINGLS_PRINGLS_20240612';
% cruise = 'NF2405_VIIRS';
% cruise = 'CHESAPEAKE_BAY_HELICOPTER_Chesapeake_Bay_2022';
% cruise = 'PVST_PRINGLS_PRINGLS_20240717';
cruise = 'PVST_PRINGLS_PRINGLS_20240813';

[fontName,projPath,imgPath] = machine_prefs();                      % <-- Set this
projPath = fullfile(projPath,'SeaBASS','JIRA_tickets',cruise);
% projPath = ...
% fullfile(projPath,'HyperPACE','field_data','DALEC',cruise,'L2','SeaBass');

if ~isfolder(projPath)
    fprintf('Bad project path: %s\n',projPath)
    return
end


% Use ls /path/path/*.sb > filelist.txt (full path filename list)
fid = fopen(fullfile(projPath,'filelist.txt')); % <-- Full path list of seabass files (JIRA)

% Es and Rrs as separate SeaBASS files
EsAndRrs = 0;                   % <--- Set to zero combined seabass files, 1 for separate Es and Rrs

%%
fileList = textscan(fid,'%s\n');
fclose(fid);

dBase = struct();

j=0;    % Allows for multiple seabass files with datetime rows ("wide" files)
for i=1:length(fileList{1})
    fp = fileList{1}{i};
    if ~contains(fp,'tgz')
        [data, sbHeader, headerArray, dataArray] = readsb(fp,'MakeStructure', true);
        fNames = fieldnames(data);
        missing = extractfield(sbHeader,'missing');

        % Need to determine the organization of the data
        if isfield(data,'wavelength')
            %%          Tall File
            disp('Found wavelength column. Single spectrum "Tall" file')
            % datenum is repeated
            dBase(i).datetime = datetime(data.datenum(1),'ConvertFrom','datenum','TimeZone','UTC');
            dBase(i).latitude = extractfield(sbHeader,'north_latitude');
            dBase(i).longitude = extractfield(sbHeader,'west_longitude');
            dBase(i).station = extractfield(sbHeader,'station');
            dBase(i).cloud = extractfield(sbHeader,'cloud_percent');
            dBase(i).wind = extractfield(sbHeader,'wind_speed');
            dBase(i).water_depth = extractfield(sbHeader,'water_depth');
            dBase(i).wave_height = extractfield(sbHeader,'wave_height');

            if isfield(sbHeader,'rho_correction')
                dBase(i).rho_correction = extractfield(sbHeader,'rho_correction');
            end
            if isfield(sbHeader,'relaz')
                dBase(i).relaz = extractfield(sbHeader,'relaz');
            end
            if isfield(sbHeader,'sza')
                dBase(i).sza = extractfield(sbHeader,'sza');
            end

            dBase(i).wavelength = data.wavelength';
            dBase(i).rrs = data.rrs';
            if isfield(data,'rrs_sd')
                dBase(i).rrs_sd = data.rrs_sd';
            elseif isfield(data,'rrs_unc')
                dBase(i).rrs_unc = data.rrs_unc';
            end

            % Es and Ed terms are assumed equivalent (Es = Ed(0+))
            if isfield(data,'es')
                dBase(i).es = data.es';
                if isfield(data,'es_sd')
                    dBase(i).es_sd = data.es_sd';
                elseif isfield(data,'es_unc')
                    dBase(i).es_unc = data.es_unc';
                end
            elseif isfield(data,'ed')
                dBase(i).es = data.ed'; % Change ed to es
                if isfield(data,'ed_sd')
                    dBase(i).es_sd = data.ed_sd';
                elseif isfield(data,'ed_unc')
                    dBase(i).es_unc = data.ed_unc';
                end
            end
            % Lref for handhelds is the plaque radiance. Check their
            % equation for Es in their notes (e.g. Es = Lref*pi*(1/0.99),
            % where 0.99 is calibrated reflectance in the relevant bands

            % Li and Lsky are the same
            if isfield(data,'li')
                dBase(i).li = data.li';
                dBase(i).li_sd = data.li_sd';
            elseif isfield(data,'lsky')
                dBase(i).li = data.lsky';
                dBase(i).li_sd = data.lsky_sd';
            end
            if isfield(data,'lw')
                dBase(i).lw = data.lw';
                dBase(i).lw_sd = data.lw_sd';
            end
            if isfield(data,'lt')
                dBase(i).lt = data.lt';
                dBase(i).lt_sd = data.lt_sd';
            end

        elseif isfield(data,'time')
            %%          Wide File
            disp('Rows organized by time. "Wide" file.')
            % This section is in early development and may not support all
            % wide submissions (e.g., unexpected fields).
            
            for n=1:length(data.time)                
                j=j+1; % Allows for multiple files with datetime rows
                dBase(j).datetime = datetime(data.datenum(n),'ConvertFrom','datenum','TimeZone','UTC');
                dBase(j).latitude = data.lat(n);
                dBase(j).longitude = data.lon(n);
                if isfield(data,'station')
                    dBase(j).station = data.station(n);
                else
                    dBase(j).station = extractfield(sbHeader,'station');
                end
                if isfield(data,'cloud')                
                    dBase(j).cloud = data.cloud(n);
                else
                    dBase(j).cloud = extractfield(sbHeader,'cloud_percent');
                end                
                if isfield(data,'wind')
                    dBase(j).wind = data.wind(n);
                else
                    dBase(j).wind = extractfield(sbHeader,'wind_speed');
                end
                dBase(j).water_depth = extractfield(sbHeader,'water_depth');
                if isfield(data,'waveht')
                    dBase(j).wave_height = data.waveht(n);
                else
                    dBase(j).wave_height = extractfield(sbHeader,'wave_height');
                end
                if isfield(data,'relaz')
                    dBase(j).relAz = data.relaz(n);
                else
                    dBase(j).relAz = extractfield(sbHeader,'relaz');
                end
                if isfield(data,'sza')
                    dBase(j).sza = data.sza(n);
                else
                    dBase(j).sza = extractfield(sbHeader,'sza');
                end                
                if isfield(data,'aot')
                    dBase(j).aot = data.aot(n);
                end

                % Spectral fields
                specFieldNames = {'rrs','ed','es','li','lt','lu','lw','rrs','lwn'};
                for s=1:length(specFieldNames)
                    sName = specFieldNames{s};

                    specCols = contains(fNames,sName,'IgnoreCase',true);
                    if any(specCols)                        
                        specUncCols = ...
                            ~cellfun(@isempty,regexp(fNames,regexptranslate('wildcard', sprintf('%s*_unc',sName)))) |...
                            ~cellfun(@isempty,regexp(fNames,regexptranslate('wildcard', sprintf('%s*_sd',sName))));

                        dataCell = struct2cell(data);
                        specCols = specCols & ~specUncCols;
                        
                        spec = horzcat(dataCell{specCols});
                        dBase(j).(sName) = spec(n,:);
                        
                        % Capture wavebands from header
                        wvFieldName = sprintf('%s_wavelength',sName);
                        dBase(j).(wvFieldName) = sbHeader.wavelength_lists.(sName);
                        
                        if any(specUncCols)
                            specUnc = horzcat(dataCell{specUncCols});
                            if strcmpi(sName,'rrs') ||strcmpi(sName,'lw') ||strcmpi(sName,'lwn')
                                % These !should! be uncertainties, not
                                % simple stds
                                fName = sprintf('%s_unc',sName);
                                dBase(j).(fName) = specUnc(n,:);
                            else
                                fName = sprintf('%s_sd',sName);
                                dBase(j).(fName) = specUnc(n,:);
                            end
                        end
                    end
                end


            end
        end
    end
end

% Nan out the missing values
fNames = fieldnames(dBase);
for col=1:length(fNames)
    dims = size(dBase(1).(fNames{col}));
    if sum(dims) == 2
        test = [dBase.(fNames{col})];
        if isa(test,'double')
            test(test == missing) = nan;
            test = num2cell(test); 
            [dBase.(fNames{col})] = test{:}; 
        end
    else
        for row=1:size(dBase,2)
            test = dBase(row).(fNames{col});
            test(test == missing) = nan;
            dBase(row).(fNames{col}) = test;
        end
    end
end

% Sort chronologically
T = struct2table(dBase);
sortedT = sortrows(T,'datetime');
dBase = table2struct(sortedT);

% Merge Es and Rrs into common dBase
if EsAndRrs
    whrRrs = find(~cellfun(@isempty,{dBase.rrs}));
    whrEs = find(~cellfun(@isempty,{dBase.es}));
    dBaseEs = dBase(whrEs);
    EsDateTime = [dBaseEs.datetime];
    for i=1:length(whrRrs)
        dBaseMerge(i) = dBase(whrRrs(i));
        dateTime1 = dBaseMerge(i).datetime;
        % Expects an exact match in timestamp
        match = find(dateTime1 == EsDateTime);
        dBaseMerge(i).es = dBaseEs(match).es;
        dBaseMerge(i).es_wavelength = dBaseEs(match).es_wavelength;
        if any(contains(fieldnames(dBase),'es_sd'))
            dBaseMerge(i).es_sd = dBaseEs(match).es_sd;
        elseif any(contains(fieldnames(dBase),'es_unc'))
            dBaseMerge(i).es_unc = dBaseEs(match).es_unc;
        end

    end
    dBase = dBaseMerge;
end

save(['dat/' cruise],'dBase')



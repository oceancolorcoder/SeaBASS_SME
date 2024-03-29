
% Aggregate all L2 HDF HyperSAS data processed in HyperCP into one
% database structure with one record per spectral ensemble. All records in
% the substructures of SASData correspond to the groups within the HDF
% files, and are related one-to-one in time and space, as reported in the
% Ancillary structure.
% Files obtained using run_get_files.bsh
%
% D. Aurin NASA/GSFC March 2024

wipe

%% Setup
cruise = 'EXPORTSNA_NASA';
% cruise = 'EXPORTSNA_Boss';
% cruise = 'viirs_2019_foster'; % Rejected for bad metadata, haha! (it's mine)


%% Run
[fontName,projPath,imgPath] = machine_prefs();
projPath = fullfile(projPath,'SeaBASS','JIRA_tickets',cruise);
outDbase = sprintf('dat/%s.mat',cruise);
disp(cruise)
SASDir = fullfile(projPath,'L2');
SASlist = dir([SASDir '/*.hdf']); % <-- Filter out STATION data? Probably.
fprintf('Number of files found: %d\n', numel(SASlist))

ensemble = 0;
for i=1:length(SASlist)
    fpf = [SASDir filesep SASlist(i).name];

    % Ignore station data for now
    if ~contains(fpf,'station','ignorecase',true)
        disp(['Loading: ', fpf])
        filename = SASlist(i).name;

        Groups = h5info(fpf).Groups; % Ancillary, Radiance, etc.
        for igrp = 1:numel(Groups)
            Datasets = h5info(fpf).Groups(igrp).Datasets;

            %% ANCILLARY
            if strcmp(Groups(igrp).Name, '/ANCILLARY')

                % Loop over datasets in Group
                for ids = 1:numel(Datasets)
                    eval(sprintf("dataset = h5read(fpf,'%s/%s/');",...
                        Groups(igrp).Name, Datasets(ids).Name))
                    dsFields = fieldnames(dataset);

                    % Filename and date/time will be common to all groups/datasets
                    % in the file. Store in Ancillary structure and
                    % eliminate elsewhere
                    if ids == 1
                        ens_start = ensemble; %Capture the output index prior to this file's first ensemble
                        eval(sprintf("dateTag = h5read(fpf, '%s/%s/').Datetag;",...
                            Groups(igrp).Name, Datasets(ids).Name))
                        eval(sprintf("timeTag = h5read(fpf, '%s/%s/').Timetag2;",...
                            Groups(igrp).Name, Datasets(ids).Name))

                        % Each dataset has multiple ensembles
                        nEns = length(dateTag);
                        for iens = 1:nEns
                            ensemble = ensemble + 1;
                            eval(sprintf("%s(ensemble).cruise = '%s';",Groups(igrp).Name(2:end), cruise))
                            eval(sprintf('%s(ensemble).filename = filename;',Groups(igrp).Name(2:end)))
                            eval(sprintf("%s(ensemble).datenum = dateTimeTag2datenum(num2str(dateTag(iens)),'%09d');", ...
                                Groups(igrp).Name(2:end), timeTag(iens) ))
                            eval(sprintf('%s(ensemble).datestring = datestr(ANCILLARY(ensemble).datenum);',...
                                Groups(igrp).Name(2:end)))
                        end
                        %                             fprintf('%s: %d ensembles\n', Groups(igrp).Name(2:end), ensemble)
                        fprintf('%s: %d ensembles\n', Groups(igrp).Name(2:end), nEns)

                        nAnc = ensemble;
                        ensemble = ens_start; % Rewind to the beginning of this file
                    end

                    for ifld = 1:numel(dsFields)
                        % Read in each non-datetime field/dataset and assign to
                        % data variable
                        if ~strcmpi(dsFields{ifld},'Datetag') && ~strcmpi(dsFields{ifld},'Timetag2')
                            eval(sprintf("data = h5read(fpf, '/%s/%s/').%s;",...
                                Groups(igrp).Name, Datasets(ids).Name, ...
                                dsFields{ifld}))
                            ens_start = ensemble; %Capture the index prior to this file's first ensemble

                            for iens = 1:numel(data)
                                ensemble = ensemble +1;

                                % For non-spectral data
                                % If the field did not already exist in the
                                % structure, it will be added with blanks
                                % before this ensemble
                                eval(sprintf('%s(ensemble).%s = data(iens);',...
                                    Groups(igrp).Name(2:end),dsFields{ifld}))
                            end

                            % Rewind to the beginning of this file for the next dataset field
                            ensemble = ens_start;
                        end
                    end
                end
            end

            %% DERIVED_PRODUCTS
            if strcmp(Groups(igrp).Name, '/DERIVED_PRODUCTS')
                for ids = 1:numel(Datasets)
                    eval(sprintf("dataset = h5read(fpf,'%s/%s/');",...
                        Groups(igrp).Name, Datasets(ids).Name))

                    % Strip off Date and Time tags
                    if isfield(dataset,'Datetag')
                        dataset = rmfield(dataset,'Datetag');
                        dataset = rmfield(dataset,'Timetag2');
                    end
                    dsFields = fieldnames(dataset);
                    nField = numel(dsFields);
                    dsFields = fieldnames(dataset);

                    % For scalar data (arbitrarily set to <4 fields in the dataset, after removing Datetag Timetag)
                    if numel(dsFields) < 4
                        for ifld = 1:nField
                            % Read in each non-datetime field/dataset and assign to
                            % data variable

                            eval(sprintf("data = h5read(fpf, '/%s/%s/').%s;",...
                                Groups(igrp).Name, Datasets(ids).Name, ...
                                dsFields{ifld}))
                            ens_start = ensemble; %Capture the index prior to this file's first ensemble

                            for iens = 1:numel(data)
                                ensemble = ensemble +1;

                                % For non-spectral data
                                eval(sprintf('%s(ensemble).%s = data(iens);',...
                                    Groups(igrp).Name(2:end),dsFields{ifld}))
                            end

                            if ids==1 && ifld==1
                                fprintf('%s: %d ensembles\n', Groups(igrp).Name(2:end), ensemble)
                                nDP = ensemble;
                                if nDP ~=nAnc
                                    disp('************************ Mis-match warning. Error! *****************')
                                end
                            end

                            % Rewind to the beginning of this file for the next
                            % field
                            if not (igrp == numel(Groups) && ids == numel(Datasets) ...
                                    && ifld == nFields && cruisei == numel(cruises))
                                ensemble = ens_start;
                            end
                        end
                    else
                        % Vector/spectral data

                        nEns = eval(sprintf('numel(dataset.%s)',dsFields{1}));
                        % Initialize an array
                        wavedata = NaN(nEns, numel(dsFields));
                        wavelength = NaN(1, numel(dsFields));
                        for wl=1:nField
                            % wavelength format here is xDDD = 3 digits D
                            k = strfind(dsFields{wl}, 'x');
                            digits = dsFields{wl}(k(1)+1:end);
                            wavelength(wl) = str2double(digits);
                            wavedata(:,wl) = dataset.(dsFields{wl});
                        end

                        ens_start = ensemble; %Capture the index prior to this file's first ensemble
                        for iens = 1:numel(data)
                            ensemble = ensemble +1;

                            % For spectral data
                            eval(sprintf('%s(ensemble).%s = wavedata(iens,:);',...
                                Groups(igrp).Name(2:end),Datasets(ids).Name))
                            eval(sprintf('%s(ensemble).%s_wavelength = wavelength;',...
                                Groups(igrp).Name(2:end),Datasets(ids).Name))
                        end

                        % Rewind to the beginning of this file for the next field
                        if not (igrp == numel(Groups) && ids == numel(Datasets)...
                                && ifld == nFields && cruisei == numel(cruises))
                            ensemble = ens_start;
                        end
                    end
                end
            end


            %% (Ir)Radiances/Reflectances
            if strcmp(Groups(igrp).Name, '/IRRADIANCE') || strcmp(Groups(igrp).Name, '/RADIANCE') || ...
                    strcmp(Groups(igrp).Name, '/REFLECTANCE')
                if strcmpi(Groups(igrp).Name,'/REFLECTANCE')
                    info = h5info(fpf,Groups(igrp).Name);
                    if contains([info.Attributes.Name],'GLINT_CORR')
                        glint = h5readatt(fpf,Groups(igrp).Name,'GLINT_CORR');
                    else
                        glint = 'None';
                    end
                    if contains([info.Attributes.Name],'NIR_RESID_CORR')
                        nir = h5readatt(fpf,Groups(igrp).Name,'NIR_RESID_CORR');
                    else
                        nir = 'None';
                    end
                end

                for ids = 1:numel(Datasets)
                    eval(sprintf("dataset = h5read(fpf,'%s/%s/');",...
                        Groups(igrp).Name, Datasets(ids).Name))
                   

                    % Strip off Date and Time tags
                    if isfield(dataset,'Datetag')
                        dataset = rmfield(dataset,'Datetag');
                        dataset = rmfield(dataset,'Timetag2');
                    end
                    dsFields = fieldnames(dataset);
                    nField = numel(dsFields);
                    dsFields = fieldnames(dataset);

                    nEns = eval(sprintf('numel(dataset.%s)',dsFields{1}));

                    for wl=length(dsFields):-1:1
                        k = strfind(dsFields{wl}, 'x');
                        if numel(k) == 1
                            % wavelength format here is xDDD = 3 digits D
                            k = strfind(dsFields{wl}, 'x');
                            digits = dsFields{wl}(k(1)+1:end);
                            wavelength(wl) = str2double(digits);
                        elseif numel(k) == 2
                            % wavelength format is xDDDDx2EP = 4 digits D and one decimal of P
                            k = strfind(dsFields{wl}, 'x');
                            digits = dsFields{wl}(k(1)+1:k(2)-1);
                            decimal = dsFields{wl}(k(2)+3);
                            digits(end) = decimal;
                            wavelength(wl) = str2double(digits)/10;
                        end
                        wavedata(:,wl) = dataset.(dsFields{wl});
                    end

                    ens_start = ensemble; %Capture the index prior to this file's first ensemble
                    for iens = 1:numel(data)
                        ensemble = ensemble +1;

                        % For spectral data
                        eval(sprintf('%s(ensemble).%s = wavedata(iens,:);',...
                            Groups(igrp).Name(2:end),Datasets(ids).Name))
                        if isempty(strfind(Datasets(ids).Name,'_unc')) && isempty(strfind(Datasets(ids).Name,'_uncorr')) ...
                                && isempty(strfind(Datasets(ids).Name,'_sd'))...
                                && isempty(strfind(Datasets(ids).Name,'_N')) && isempty(strfind(Datasets(ids).Name,'nir_')) ...
                                && isempty(strfind(Datasets(ids).Name,'rho_'))
                            % && isempty(strfind(Datasets(ids).Name,'_median')) && isempty(strfind(Datasets(ids).Name,'_sd'))...

                            eval(sprintf('%s(ensemble).%s_wavelength = wavelength;',...
                                Groups(igrp).Name(2:end),Datasets(ids).Name))
                        end
                    end
                    clear wavelength wavedata

                    if ids==1
                        fprintf('%s: %d ensembles\n', Groups(igrp).Name(2:end), nEns)
                    end

                    % Rewind to the beginning of this file for the next field
                    if not (igrp == numel(Groups) && ids == numel(Datasets))
                        ensemble = ens_start;
                    end
                end

            end
        end
    end
    fprintf(' Cruise running Total: %d ensembles\n', size(ANCILLARY,2))
end

SASData.Ancillary = ANCILLARY;
if exist('DERIVED_PRODUCTS','var')
    SASData.Derived_Products = DERIVED_PRODUCTS;
end
SASData.Radiance = RADIANCE;
SASData.Irradiance = IRRADIANCE;
SASData.Reflectance = REFLECTANCE;
clear ANCILLARY DERIVED_PRODUCTS RADIANCE IRRADIANCE REFLECTANCE

fprintf(' Database running Total: %d ensembles\n', size(SASData.Ancillary,2))

% Sort by date
[~,idx]=sort([SASData.Ancillary.datenum]);

SASData.Ancillary = SASData.Ancillary(idx);
SASData.Reflectance = SASData.Reflectance(idx);
SASData.Irradiance = SASData.Irradiance(idx);
SASData.Radiance = SASData.Radiance(idx);
if exist('DERIVED_PRODUCTS','var')
    SASData.Derived_Products = SASData.Derived_Products(idx);
end

save(outDbase,'SASData')

%%
% function mergedGroup = catGroup(oldGroup, newGroup)
% % This is mental...
% % Concatonate the new group to the old group making sure they have the
% % same fields in the same order.
% nEnsOld = numel(oldGroup);
% nEnsNew = numel(newGroup);
% oldAncFields = fieldnames(oldGroup);
% newAncFields = fieldnames(newGroup);
% for iflds = 1:numel(oldAncFields)
%     if ~(ismember(oldAncFields{iflds},newAncFields))
%         temp = num2cell(NaN(nEnsNew,1));
%         eval(sprintf('[newGroup(:).%s] = deal(temp{:});', oldAncFields{iflds}))
%     end
% end
% for iflds = 1:numel(newAncFields)
%     if ~(ismember(newAncFields{iflds},oldAncFields))
%         % This is mental...
%         temp = num2cell(NaN(nEnsOld,1));
%         eval(sprintf('[oldGroup(:).%s] = deal(temp{:});', newAncFields{iflds}))
%     end
% end
% oldAncFields = fieldnames(oldGroup);
% newGroup = orderfields(newGroup,oldAncFields);
% mergedGroup = [oldGroup, newGroup];
% end

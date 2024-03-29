% Using readsb.m (in MATLAB_embedded, updated No 2023), pull in all the
% files obtained using run_get_files.bsh.
%
% Output:
%   dBase structure with relevant fields (e.g., datetime, Rrs, etc.)
%
% D. Aurin NASA/GSFC March 2024

wipe
% cruise = 'viirs_2019_foster'; % Use make_database_hypercp.m
% cruise = 'RSWQ_2023'; % single spectrum per file
% cruise = 'JackBlanton'; % Rivero-Calle; returned to PI
cruise = 'Belgium_2021'; % Twardowski; Not AWR. Returned to Mike.

[fontName,projPath,imgPath] = machine_prefs();
projPath = fullfile(projPath,'SeaBASS','JIRA_tickets',cruise);

fid = fopen(fullfile(projPath,'filelist.txt'));
fileList = textscan(fid,'%s\n');
fclose(fid);

dBase = struct();

for i=1:length(fileList{1})
    file = fileList{1}{i};
    if ~contains(file,'tgz')

        % fp = fullfile(projPath,file);
        fp = file;
        [data, sbHeader, headerArray, dataArray] = readsb(fp,'MakeStructure', true);
        fNames = fieldnames(data);

        % Test for multiple spectra per file
        % if length(data.datenum)

        % Need to determine the organization of the data
        if sum(strcmpi(fNames,'wavelength')) ~= 0
            disp('Found wavelength column. Single spectrum file')            
            % datenum is repeated
            dBase(i).datetime = datetime(data.datenum(1),'ConvertFrom','datenum','TimeZone','UTC');
            dBase(i).latitude = extractfield(sbHeader,'north_latitude');
            dBase(i).longitude = extractfield(sbHeader,'west_longitude');
            dBase(i).station = extractfield(sbHeader,'station');
            dBase(i).cloud = extractfield(sbHeader,'cloud_percent');
            dBase(i).wind = extractfield(sbHeader,'wind_speed');
            dBase(i).water_depth = extractfield(sbHeader,'water_depth');
            dBase(i).rrs = data.rrs';
            dBase(i).rrs_sd = data.rrs_sd';
            dBase(i).es = data.ed';
            dBase(i).es_sd = data.ed_sd';
            dBase(i).wavelength = data.wavelength';

            missing = extractfield(sbHeader,'missing');

        else
            if sum(strcmpi(fNames,'time')) ~= 0
                if length(data.time) > 1
                    disp('Multiple timestamps found. Multi-spectrum file.')
                end
            end
        end
    end
end

% Nan out the nans
fNames = fieldnames(dBase);
for col=1:length(fNames)
    dims = size(dBase(1).(fNames{col}));
    if sum(dims) == 2
        test = [dBase.(fNames{col})];
        if isa(test,'double')
            test(test == missing) = nan;
            test = num2cell(test); % Ack!
            [dBase.(fNames{col})] = test{:}; % Ack!
        end
    else
        for row=1:size(dBase,2)
            test = dBase(row).(fNames{col});
            test(test == missing) = nan;
            dBase(row).(fNames{col}) = test;
        end
    end
end

save(['dat/' cruise],'dBase')
    


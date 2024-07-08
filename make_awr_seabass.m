% Using readsb.m (in MATLAB_embedded, updated No 2023), pull in all the
% files obtained using run_get_files.bsh.
%
% Output:
%   dBase structure with relevant fields (e.g., datetime, Rrs, etc.)
%
% D. Aurin NASA/GSFC March 2024

wipe
% cruise = 'viirs_2019_foster'; 
% cruise = 'RSWQ_2023'; % single spectrum per file
% cruise = 'JackBlanton'; % Rivero-Calle; returned to PI
% cruise = 'Belgium_2021'; % Twardowski; Not AWR. Returned to Mike.
% cruise = 'UMCES_Missouri_Reservoirs'; % SBA Lorena Silva/Greg Silsbe
% cruise = 'BIOSCAPE_COASTAL_CARBON_Walker_Bay'; % Kyle Turner/Maria Tzortziou, BIOSCAPE (S. Africa)
% cruise = 'Brewin_Superyacht_Science_2018';
cruise = 'ArcticCC_Norton_Sound_2022';

[fontName,projPath,imgPath] = machine_prefs();
projPath = fullfile(projPath,'SeaBASS','JIRA_tickets',cruise);

fid = fopen(fullfile(projPath,'filelist.txt'));
fileList = textscan(fid,'%s\n');
fclose(fid);

dBase = struct();

j=0;
for i=1:length(fileList{1})
    fp = fileList{1}{i};
    if ~contains(fp,'tgz')
        [data, sbHeader, headerArray, dataArray] = readsb(fp,'MakeStructure', true);
        fNames = fieldnames(data);

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
            dBase(i).wave_height = extractfield(sbHeader,'wave_height');
            dBase(i).wavelength = data.wavelength';
            dBase(i).rrs = data.rrs';
            dBase(i).rrs_sd = data.rrs_sd';
            % Es and Ed terms are assumed equivalent (Es = Ed(0+))
            if sum(contains(fieldnames(data),'es')) > 0
                dBase(i).es = data.es';
                dBase(i).es_sd = data.es_sd';
            elseif sum(contains(fieldnames(data),'ed')) > 0
                dBase(i).es = data.ed';
                dBase(i).es_sd = data.ed_sd';
            end
            % Lref for handhelds is the plaque radiance. Check their
            % equation for Es in their notes (e.g. Es = Lref*pi*(1/0.99),
            % where 0.99 is calibrated reflectance in the relevant bands

            % Li and Lsky are the same
            if sum(contains(fieldnames(data),'li')) > 0
                dBase(i).li = data.li';
                dBase(i).li_sd = data.li_sd';
            elseif sum(contains(fieldnames(data),'lsky')) > 0
                dBase(i).li = data.lsky';
                dBase(i).li_sd = data.lsky_sd';
            end
            if sum(contains(fieldnames(data),'lw')) > 0
                dBase(i).lw = data.lw';
                dBase(i).lw_sd = data.lw_sd';
            end
            if sum(contains(fieldnames(data),'lt')) > 0
                dBase(i).lt = data.lt';
                dBase(i).lt_sd = data.lt_sd';
            end
            

            missing = extractfield(sbHeader,'missing');

        else
            if sum(strcmpi(fNames,'time')) ~= 0                
                disp('Rows organized by time')
                if j==0
                    dBase(1).datetime = datetime(now,'ConvertFrom','datenum','TimeZone','UTC');
                end
                % For HyperInSPACE output, this will require matching Es to
                % Rrs from seperate files
                for n=1:length(data.time)
                    % Check for existing matching time in dBase
                    newDateTime = datetime(data.datenum(n),'ConvertFrom','datenum','TimeZone','UTC');                    
                    [x,y] = find_nearest(newDateTime,[dBase.datetime]);
                    if x ~= newDateTime
                        j=j+1;
                        dBase(j).datetime = datetime(data.datenum(n),'ConvertFrom','datenum','TimeZone','UTC');
                        dBase(j).latitude = datetime(data.latitude(n));
                        dBase(j).latitude = datetime(data.latitude(n));
                        dBase(i).station = extractfield(sbHeader,'station');
                        % dBase(i).cloud = extractfield(sbHeader,'cloud_percent');
                        dBase(i).cloud = data.cloud(n);
                        % dBase(i).wind = extractfield(sbHeader,'wind_speed');
                        dBase(i).wind = data.wind(n);
                        dBase(i).water_depth = extractfield(sbHeader,'water_depth');
                        dBase(i).wave_height = extractfield(sbHeader,'wave_height');                        
                        dBase(i).relAz = data.relAz(n);
                        dBase(i).sza = data.sza(n);
                        dBase(i).aot = data.aot(n);

                        

                    else                        
                        dBase(y).dataXXX = XXX;
                    end

                    

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

% Sort chronologically
T = struct2table(dBase);
sortedT = sortrows(T,'datetime');
dBase = table2struct(sortedT);

save(['dat/' cruise],'dBase')
    


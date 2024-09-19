function flags = writeFlags(ancillary, thresholds, AWR)
% flags is the structure with individual flags
% flag is the logical vector for the prescribed flags
% 0, 1, or 2 for reject, seabass-only, validation


%% Set flags
if ancillary.validation
    flags.Cloud = (ancillary.cloud > thresholds.cloud(1)) & (ancillary.cloud < thresholds.cloud(2)); %Partly cloudy
else
    flags.Cloud = false(1,length(ancillary.cloud));
end
flags.Wind = ancillary.wind > thresholds.wind;
flags.SZA = ancillary.sza < thresholds.sza(1) | ancillary.sza > thresholds.sza(2);

if ~ancillary.SBA
    flags.RelAz = abs(ancillary.relAz)<thresholds.relAz(1) | abs(ancillary.relAz)>thresholds.relAz(2);
end
flags.QWIP = ancillary.qwip > thresholds.qwip;
flags.QA = ancillary.qa < thresholds.qa;
waveRange = find(AWR.wave >= thresholds.negRrs(1) & AWR.wave <= thresholds.negRrs(2));
flags.negRrs = any(AWR.rrs(:,waveRange) < 0.0,2)';

flags.Manual = 0*flags.negRrs;

save(sprintf('dat/%s_flags.mat',ancillary.cruise),"AWR","ancillary","flags")


%% Apply flags
if ~ancillary.SBA
    flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.RelAz] | [flags.QWIP] ...
        | [flags.negRrs];
else
    % No RelAz flag
    flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.QWIP] ...
        | [flags.negRrs];
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
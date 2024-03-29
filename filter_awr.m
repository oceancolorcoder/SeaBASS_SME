% Use the flags and ancillary data from review_awr_[type].m to rank the spectra
%
% Rank as   1) Validation candidate
%           2) Validation candidate with changes
%           3) SeaBASS only
%           4) Fail
%   These designations are subjective, but technically for 1 it would need
%   to be cloud free which would have to come from ancillary data. The
%   cloud flag is set for identifying anything between 20%-80% cloud (fully
%   cloudy should not impact glint correction, though not useful for validation).
%       If validation is set to only cloud=0%, nearly all data would be
%       lost.
%
%   Validation candidates must pass all flags.
%   Identify any spectra unworthy of SeaBASS:
%       QWIP flag + Visual inspection
%
% D. Aurin NASA/GSFC March 2024

%% Setup
wipe
fontName = machine_prefs;
cruise = 'EXPORTSNA_NASA';
% cruise = 'EXPORTSNA_Boss';

%%
load(sprintf('dat/%s_flags.mat',cruise))

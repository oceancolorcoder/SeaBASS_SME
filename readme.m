
%% run_get_files.bsh
%   Inputs: 
%       [filelistname].lis from Jira
%       provide a unique cruisename (does not have to match SeaBASS)
%   
%   Downloads all the SeaBASS files and moves them to a cruise directory
%   under Projects/SeaBASS/Jira_tickets

%% make_database_hypercp.m
%   Inputs:
%       cruisename
%       L2 HDF5 files from HyperCP
%
%   Outputs:
%       dat/[cruisename].mat database containing SASData structure
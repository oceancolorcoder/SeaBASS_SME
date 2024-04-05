## run_get_files.bsh
    Inputs:
       [filelistname].lis from Jira
       provide a unique cruisename (does not have to match SeaBASS)

   Downloads all the SeaBASS files and moves them to a cruise directory
   under Projects/SeaBASS/Jira_tickets

## make_awr_seabass.m
    Inputs:
       cruisename
       SeaBASS files from .lis

    Outputs:
       dat/[cruisename].mat database containing dBase structure

## make_awr_hypercp.m
    Inputs:
       cruisename
       L2 HDF5 files from HyperCP

    Outputs:
       dat/[cruisename].mat database containing SASData structure

## review_awr_hypercp.m
    Inputs:
       dat/[cruisename].mat from make_awr_hypercp.m (SASData structure)

    Outputs:
       dat/[cruisename]_flags.mat (select variables and QA/QC flags)
       plt/[cruisename]_AllSpec.png (Rrs, Es, Li, Lw with flags shown)
       plt/[cruisename]_QWIP.png
       plt/[cruisename]_QWIP_Hist.png
       plt/[cruisename]_timeline.png (time series of select variables used
       in QAQC)

 Add AVW and QWIP if not already included in Derived Products
 Display spectral (ir)radiance and reflectance with ancillary data
 dashboard (timeline)\

 Flag for QWIP,RelAz,Wind,SZA, negative Rrs, and visual inspection\
 Thresholds\
 For .env.all:\
 negRrs = [380 680]; % Spectral range of negatives to eliminate from all sets\
 tRelAz = [87 138]; % M99, Z17, IOCCG\
 tSZA = [18 62]; % e.g. 20: Zhang 2017, depends on wind, e.g. 60:Brewin 2016\
 tWind = 10; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,\
 tQWIP = 0.2; % Dierssen et al. 2022\
 tQA = 0.2; % This is more experimental. Monitor it, but don't filter it.\
 tCloud = [20 80]; % Clear and fully overcast should be okay. 20% - 80% are likely to be the worst. This is experimental.\
 cloudIndexes = [0.05 0.3]; % From Ruddick et al. 2006 based on M99 models, where <0.05 is clear, >0.3 is fully overcast\

% Thresholds for validation\
tRelAz = [89 136]; % M99, Z17, IOCCG\
tWind = 7; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,\
tQWIP = 0.17; % Experimental\

## review_awr_seabass.m
    Inputs:
       dat/[cruisename].mat from make_awr_seabass.m (dBase structure)

   Otherwise the same as review_awr_hypercp.m

## batch_awr2env.py
    Inputs:
        List of .sb files to process

    Outputs:
        .env.all files
Runs awr2env.py

## awr2env.py
    Inputs:
        SeaBASS files identified in batch_awr2env.py (Es, Rrs)

    Outputs:
        *.env.all file for NOMAD
        *.env file for VALIDATION
Matches spectra to within 10s of flags developed in review_awr_TYPE.m and only includes spectra that passed all flags.



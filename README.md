## run_get_files.bsh
    Inputs:
       [filelistname].lis from Jira
       provide a unique cruisename (does not have to match SeaBASS)

   Downloads all the SeaBASS files and moves them to a cruise directory
   under Projects/SeaBASS/Jira_tickets

## make_awr_seabass.m
    Inputs:
       cruisename
       SeaBASS files from .list

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
 dashboard (timeline)

 Flag for QWIP,RelAz,Wind,SZA, and visual inspection
 Thresholds
 tRelAz = [88 137];  M99, Z17, IOCCG
 tWind = 10;   6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,
 tSZA = [18 62];  e.g. 20: Zhang 2017, depends on wind, e.g. 60:Brewin 2016
 tQWIP = 0.2;  Dierssen et al. 2022
 tQA = 0.2;  This is more experimental. Monitor it, but don't filter it.
 tCloud = [20 80];  Clear and fully overcast should be okay. 20 - 80 are likely to be the worst. This is experimental.
 cloudIndexes = [0.05 0.3];  From Ruddick et al. 2006 based on M99 models, where <0.05 is clear, >0.3 is fully overcast

## review_awr_seabass.m
    Inputs:
       dat/[cruisename].mat from make_awr_seabass.m (dBase structure)

   Otherwise the same as review_awr_hypercp.m

## awr2env.py
    Inputs:
        SeaBASS files identified in batch_awr2env.py

    Outputs:
        *.env.all file

A work in progress. Currently accommodates hyperspectral Es, Rrs, but overwrites the same [cruise].awr.[pi].env.all file with every file in the input list. Also, need to introduce timestamp screening to only output data for validation, non-validation, etc.

## batch_awr2env.py
    Inputs:
        List of .sb files to process

    Outputs:
        .env.all files



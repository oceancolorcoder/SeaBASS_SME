#!/bin/bash
projDir='/Users/daurin/Projects/SeaBASS/JIRA_tickets'
cruise='UMCES_Missouri_Reservoirs'
for file in ${projDir}/${cruise}/*.*; do
    if [[ $file == *.tgz ]];then
        echo 'Renaming/Untarzipping ' $file
        newFilename=${file/.sb/}
        tar -xvf $newFilename -C ${projDir}/${cruise}
    fi
done
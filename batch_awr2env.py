'''
SeaBASS_SME tool for converting .sb files to .env files
using awr2env.py
2024-05-05 D. Aurin NASA/GSFC
'''

import sys
import os
from pathlib import Path
import glob
import shutil
import awr2env_wide
import awr2env_tall
import awr2env_tall_thomas

# cruise='EXPORTS_EXPORTSNA_Boss_AOP_pySAS_R2'
# cruise='EXPORTS_EXPORTSNA_Mannino_AOP_pySAS_R1'
# cruise='EXPORTS_EXPORTSNP_Mannino_AOP_HyperSAS_R0'
# cruise='UMCES_Missouri_Reservoirs'
# cruise='Brewin_Superyacht_Science_2019-2020'
# cruise='Brewin_Superyacht_Science_2018'
# cruise='BIOSCAPE_COASTAL_CARBON_Walker_Bay'
# cruise='BIOSCAPE_COASTAL_CARBON_St_Helena_Bay_2023'
# cruise='VIIRS_VALIDATION_viirs_2021_gunter'
# cruise='VIIRS_VALIDATION_viirs_2022_sette'
# cruise='VIIRS_VALIDATION_viirs_2023_shimada'
# cruise='KORUS_KR_2016_RV_Onnuri_HyperSAS'
# cruise='NORTHERN_INDIAN_OCEAN_EKAMSAT-EKAMSAT-2024-Bay-of-Bengal'
# cruise='PVST_PRINGLS_PRINGLS_20240417'
# cruise='PVST_PRINGLS_PRINGLS_20240513'
# cruise='PVST_PRINGLS_PRINGLS_20240612'
# cruise='CHESAPEAKE_BAY_HELICOPTER_Chesapeake_Bay_2022'
# cruise='PVST_PRINGLS_PRINGLS_20240717'
# cruise='PVST_PRINGLS_PRINGLS_20240813'
# cruise='PVST_PRINGLS_PRINGLS_20241003'
cruise='PVST_PRINGLS_PRINGLS_20240911'

tall_thomas = False

metadata={}
metadata['dataType']='AOP'

# metadata['instrument']='HyperCP'
# metadata['instrument']='TriOS'
metadata['instrument']='SVC'

# metadata['subInstrument']='pySAS'
# metadata['subInstrument']='SolarTracker'
# metadata['subInstrument']='SBA'
metadata['subInstrument']='SVC'
# metadata['subInstrument']='GER'


if sys.platform == 'darwin':
    basePath = '/Users/daurin/Projects/SeaBASS/JIRA_tickets'
else:
    basePath = '/accounts/daurin/Projects/SeaBASS/JIRA_tickets'

# inPath = os.path.join(basePath,cruise,'test')
inPath = os.path.join(basePath,cruise)

for all_flags in [True,False]:
    # for type in ['Es','Rrs']:
    for rad_type in ['Rrs']:
        if metadata['instrument'].lower() == 'hypercp':
            # non-STATIONs from HyperCP
            fileList = glob.glob(os.path.join(inPath,f'*{rad_type}_[!STATION]*.sb'))
        elif metadata['instrument'].lower() == 'svc':
            fileList = glob.glob(os.path.join(inPath,'*.sb'))
        else:
            fileList = glob.glob(os.path.join(inPath,f'*{rad_type}_*.sb'))

        if not fileList:
            print('File list is empty. Check path.')
            print(f'{inPath}')
            break

       # call(['python','awr2env.py','--seabass_file'] + fileList)
        dict_args={}
        dict_args['seabass_file'] = fileList
        dict_args['metadata'] = metadata
        if all_flags:
            dict_args['flag_file'] = f'./dat/{cruise}_all_flags.csv'
            dict_args['all'] = True
        else:
            dict_args['flag_file'] = f'./dat/{cruise}_flags.csv'
            dict_args['all'] = False

        # Convert cruise to .env or .env.all
        if metadata['instrument'].lower() == 'hypercp':
            # For files with no wavelength column
            fileout_env = awr2env_wide.main(dict_args)
        else:
            if tall_thomas:
                fileout_env = awr2env_tall_thomas.main(dict_args)
            else:
                fileout_env = awr2env_tall.main(dict_args)
        source = Path('./',fileout_env)
        # dest = Path('./dat/',fileout_env)
        dest = Path(inPath,'ENV')
        if not dest.exists():
            Path.mkdir(dest)
        if source.exists():
            shutil.move(os.path.join('.',fileout_env),os.path.join(dest,fileout_env))


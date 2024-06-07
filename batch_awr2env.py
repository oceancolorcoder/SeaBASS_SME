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
import awr2env_hypercp
import awr2env_seabass

# cruise='EXPORTSNA_Boss'
# cruise='EXPORTSNA_NASA'
cruise='UMCES_Missouri_Reservoirs'

metadata={}
metadata['dataType']='AOP'
# metadata['instrument']='HyperCP'
metadata['instrument']='TriOS'
# metadata['subInstrument']='pySAS'
metadata['subInstrument']='SBA'


if sys.platform == 'darwin':
    basePath = '/Users/daurin/Projects/SeaBASS/JIRA_tickets'
else:
    basePath = '/accounts/daurin/Projects/SeaBASS/JIRA_tickets'

# inPath = os.path.join(basePath,cruise,'test')
inPath = os.path.join(basePath,cruise)

for all in [True,False]:
    # for type in ['Es','Rrs']:
    for type in ['Rrs']:
        if metadata['instrument'].lower() == 'hypercp':
            # STATIONs from HyperCP
            fileList = glob.glob(os.path.join(inPath,f'*{type}_[!STATION]*.sb'))
        else:
            fileList = glob.glob(os.path.join(inPath,f'*{type}_*.sb'))

        # call(['python','awr2env.py','--seabass_file'] + fileList)
        dict_args={}
        dict_args['seabass_file'] = fileList
        dict_args['metadata'] = metadata
        if all:
            dict_args['flag_file'] = f'./dat/{cruise}_all_flags.csv'
            dict_args['all'] = True
        else:
            dict_args['flag_file'] = f'./dat/{cruise}_flags.csv'
            dict_args['all'] = False

        # Convert cruise to .env or .env.all
        if metadata['instrument'].lower() == 'hypercp':
            # For files with no wavelength column
            fileout_env = awr2env_hypercp.main(dict_args)
        else:
            fileout_env = awr2env_seabass.main(dict_args)
        source = Path('./',fileout_env)
        dest = Path('./dat/',fileout_env)
        if source.exists():
            shutil.move(fileout_env,dest)




'''
SeaBASS_SME tool for converting .sb files to .env files
using awr2env.py
2024-03-29 D. Aurin NASA/GSFC
'''

import sys
import os
import glob
from subprocess import call
# import awr2env

cruise='EXPORTSNA_NASA'
if sys.platform == 'darwin':
    basePath = '/Users/daurin/Projects/SeaBASS/JIRA_tickets'
else:
    basePath = '/accounts/daurin/Projects/SeaBASS/JIRA_tickets'
inPath = os.path.join(basePath,cruise,'test')

type='Es'
fileList = glob.glob(os.path.join(inPath,f'*{type}*.sb'))

call(['python','awr2env.py','--seabass_file'] + fileList)
# call(['python','awr2env.py'])
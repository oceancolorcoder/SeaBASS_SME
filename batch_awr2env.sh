#!/bin/bash
'''
SeaBASS_SME tool for converting .sb files to .env files
using awr2env.py
2024-03-29 D. Aurin NASA/GSFC
'''

cruise='EXPORTSNA_NASA'
basePath = '/Users/daurin/Projects/SeaBASS/JIRA_tickets/test'
for file in ${basePath}/*.*; do
    python awr2env.py file
done

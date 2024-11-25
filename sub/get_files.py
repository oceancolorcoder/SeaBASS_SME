#!/usr/bin/env python3
# coding: utf-8

'''
written by jpscott6 on 2018/06/18 (joel.scott@nasa.gov)
'''

def main():

    import argparse
    import os
    import re

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
      This program downloads a list of files from a given input text file.

      REQUIRED inputs:
          1) --ifile=        a ASCII text file containing a list of file URLs

      Outputs:
          1) downloaded files to the current working directory

      Example usage call:
         get_files.py --ifile=[file name].txt

      Caveats:
        * Compatibility: This script was developed for Python 3

      License:
        /*=====================================================================*/
                         NASA Goddard Space Flight Center (GSFC) 
                 Software distribution policy for Public Domain Software

         The fd_matchup.py code is in the public domain, available without fee for 
         educational, research, non-commercial and commercial purposes. Users may 
         distribute this code to third parties provided that this statement appears
         on all copies and that no charge is made for such copies.

         NASA GSFC MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THE SOFTWARE
         FOR ANY PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED
         WARRANTY. NEITHER NASA GSFC NOR THE U.S. GOVERNMENT SHALL BE LIABLE FOR
         ANY DAMAGE SUFFERED BY THE USER OF THIS SOFTWARE.
        /*=====================================================================*/
      ''',add_help=True)

    parser.add_argument('--ifile', nargs=1, required=True, type=str, help='''\
      REQUIRED: input ASCII text file containing a list of file URLs to download
      ''')

    args=parser.parse_args()
    if not args.ifile:
        parser.error("please specify a --ifile")
    else:
        dict_args=vars(args)

    # input verification
    if not dict_args['ifile'][0] or not os.path.isfile(dict_args['ifile'][0]):
        parser.error('invalid --ifile specified, does ' + dict_args['ifile'][0] + ' exist?')
    
    fin = open(dict_args['ifile'][0], 'r')
    urls = fin.readlines()
    fin.close()

    urls = [re.sub("[\r\n]+",'',url).strip() for url in urls]
    
    for url in urls:
        if not url:
            continue
        else:
            print('Downloading:',url)
            fname = download_file(url, './')
            print('Successfully fetched:',fname)

#==========================================================================================================================================

def download_file(url, out_dir):
    '''
    download_file downloads a file
    given URL and out_dir strings
    syntax fname_local = download_file(url, out_dir)
    '''
    from pathlib import Path
    import requests
    import os

    out_path = Path(out_dir).expanduser() #even if already a Path object, this convert should work
    local_filename = out_path / url.split('/')[-1]

    try:
        r = requests.get(url, stream=True)
        r.raise_for_status()

    except requests.exceptions.HTTPError as err:
        print('Error in download_file:',err)
        return -1

    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024): 
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)

    if not local_filename.exists():
        print('Error in download_file: local copy of file not found')
        return -2

    return str(local_filename)

#==========================================================================================================================================

if __name__ == "__main__": main()

#!/usr/bin/env python3
# coding: utf-8

"""
A Perl script to create an env file(s) for data file(s) listed containing bottle data
adapted from PERL script (/disk03/PERL/scripts/ENV/bottle2env.pl)
by J.Scott on 2018/11/02 (joel.scott@nasa.gov)
Adapted for Above Water Radiometry (AWR) by D. Aurin 2024-03-29
"""

def main():

    import argparse
    from pathlib import Path
    from datetime import datetime
    from math import isnan
    import re
    from copy import copy
    from collections import OrderedDict
    from SB_support import readSB#, is_number

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
      This program creates an env.all file(s) for a given SeaBASS data file(s)
      containing Above Water Radiometry (AWR) and outputting a standard set of
      headers and fields common to the ENV file format.

      ENV.all files are interim files that are then culled into ENV files which are
      used to define and load validation matchup targets and data into the SeaBASS
      matchup work-flow.

      REQUIRED inputs:
          1) --seabass_file=    a valid SeaBASS file(s) with latitude, longitude, and date-time, and one or more of:
                                Rrs, Es.

      Outputs:
          1) an ENV.all file with a reduced set of headers and only: /fields=year,month,day,hour,minute,second,lat,lon,depth,
             and one or more of: Rrs, Es
             in the current working directory with a filename of the form: <cruise_name>.awr.<PI_name>.env.all

      Example usage call:
         awr2env.py --seabass_file filename1.sb [filename2.sb filename3.sb ...]
      ''',add_help=True)

    parser.add_argument('--seabass_file', nargs='+', type=argparse.FileType('r'), required=True, help='''\
      REQUIRED: input SeaBASS file(s)
      Must be a valid SeaBASS file containing Rrs, or Es, and latitude, longitude, and date-time.
      ''')

    args=parser.parse_args()
    dict_args=vars(args)

    fields_req = ['rrs', 'es']

    fields_opt = ['bincount']

    # fields_dep = ['depth', 'pressure']

    missing = '-9999'

    #loop over input SB files
    for fnamein in dict_args['seabass_file']:

        filein_sb = Path(fnamein.name)

        # read and verify SeaBASS file and required fields
        if filein_sb.exists():

            ds = readSB(filename=str(filein_sb), \
                        mask_missing=True, \
                        mask_above_detection_limit=True, \
                        mask_below_detection_limit=True, \
                        no_warn=True)

        else:

            parser.error('ERROR: invalid --seabass_file specified; does ' + filein_sb.name + ' exist?')

        #check for valid cruise name
        if not 'cruise' in ds.headers:

            parser.error('ERROR: cruise name not found in ' + filein_sb.name)

        #check for valid fields
        if not 'fields' in ds.headers:

            parser.error('ERROR: fields not found in ' + filein_sb.name)

        #check for valid dtimes
        ds.dtime = ds.fd_datetime()

        if not ds.dtime:

            parser.error('ERROR: date-time not parsable in ' + filein_sb.name)

        # header_depth = False
        # depth_field  = ''
        fields_fou = []

        for var in ds.variables.keys():

            # #check for depth fields in data matrix
            # for field in fields_dep:

            #     ma = re.search('^' + field + '$', var)

            #     if ma:

            #         depth = copy(missing)
            #         depth_field = field

            #check for required fields
            for field in fields_req:

                # ma = re.search('^' + field + '$', var)
                ma = re.search('^' + field + '[0-9][0-9][0-9]', var)

                if ma:

                    # fields_fou.append(field)
                    fields_fou.append(var)

                #check for error fields
                for esuf in ds.err_suffixes:

                    efield = field + esuf

                    ma = re.search('^' + efield + '$', var)

                    if ma:

                        fields_fou.append(efield)

            #check for optional fields
            for field in fields_opt:

                ma = re.search('^' + field + '$', var)

                if ma:

                    fields_fou.append(field)


        if not fields_fou:

            parser.error('ERROR: AWR data not found in ' + filein_sb.name)

        # if not depth:

        #     if 'measurement_depth' in ds.headers:

        #         if ds.headers['measurement_depth'] > 0.0:

        #             header_depth = True
        #             depth = ds.headers['measurement_depth']

        #         else:

        #             parser.error('ERROR: valid depth data not found in data matrix nor header in ' + filein_sb.name)

        #     else:

        #         parser.error('ERROR: depth data not found in data matrix nor header in ' + filein_sb.name)

        # define output vars
        out_dir = Path('./')

        fileout_sb = ds.headers['cruise'].lower() + '.awr.' + ds.pi.lower() + '.env.all'

        lat_lis  = []
        lon_lis  = []
        unit_out  = OrderedDict()
        data_out = OrderedDict()

        # fill output vars
        data_out['dt'] = []

        data_out['year'] = []
        unit_out['year'] = 'yyyy'

        data_out['month'] = []
        unit_out['month'] = 'mo'

        data_out['day'] = []
        unit_out['day'] = 'dd'

        data_out['hour'] = []
        unit_out['hour'] = 'hh'

        data_out['minute'] = []
        unit_out['minute'] = 'mn'

        data_out['second'] = []
        unit_out['second'] = 'ss'

        data_out['lat'] = []
        unit_out['lat'] = 'degrees'

        data_out['lon'] = []
        unit_out['lon'] = 'degrees'

        # data_out[depth_field] = []
        # unit_out[depth_field] = ds.variables[depth_field][1]

        data_out['associated_files'] = []
        unit_out['associated_files'] = 'none'

        data_out['associated_file_types'] = []
        unit_out['associated_file_types'] = 'none'

        for i in range(ds.length):

            #verify each row of lat, lon
            if isnan(float(ds.data['lat'][i])) or \
                isnan(float(ds.data['lon'][i])):

                continue

            # #verify each row of depth
            # if not header_depth:

            #     if isnan(float(ds.data[depth_field][i])):

            #         continue

            #     else:

            #         depth = ds.data[depth_field][i]

            #verify row has valid AWR data
            flag_nodat = 0

            for field in fields_fou:

                if isnan(ds.data[field][i]):

                    flag_nodat = flag_nodat + 1

            if flag_nodat == len(fields_fou):

                continue

            #fill meta-data
            data_out['dt'].append(ds.dtime[i])

            data_out['year'].append(ds.dtime[i].strftime('%Y'))
            data_out['month'].append(ds.dtime[i].strftime('%m'))
            data_out['day'].append(ds.dtime[i].strftime('%d'))
            data_out['hour'].append(ds.dtime[i].strftime('%H'))
            data_out['minute'].append(ds.dtime[i].strftime('%M'))
            data_out['second'].append(ds.dtime[i].strftime('%S'))

            lat_lis.append(ds.data['lat'][i])
            lon_lis.append(ds.data['lon'][i])

            data_out['lat'].append('{:.4f}'.format(ds.data['lat'][i]))
            data_out['lon'].append('{:.4f}'.format(ds.data['lon'][i]))

            # data_out[depth_field].append(depth)

            data_out['associated_files'].append(filein_sb.name)
            data_out['associated_file_types'].append('env')

            #fill required data into output dict
            for field in fields_fou:

                #handle field
                if not field in data_out:

                    data_out[field] = []
                    unit_out[field] = ds.variables[field][-1]

                if not isnan(ds.data[field][i]):

                    data_out[field].append(ds.data[field][i])

                else:

                    data_out[field].append('nan')

        #sort data_out by dtime
        for field in data_out:

            if not 'dt' in field:

                data_out[field] = [x for (y,x) in sorted(zip(data_out['dt'],data_out[field]), key=lambda pair: pair[0])]

        # create and fill fileout_sb
        print('Creating', fileout_sb, 'from:', filein_sb.name)

        with open(out_dir / fileout_sb, 'w') as fout:

            #output headers
            fout.write('/begin_header\n')
            fout.write('/cruise=' + ds.headers['cruise'].lower() + '\n')
            fout.write('/data_file_name=' + fileout_sb + '\n')
            fout.write('/data_type=matchup\n')

            #output date-time headers
            dt_min = min(data_out['dt'])
            dt_max = max(data_out['dt'])
            fout.write('/start_date=' + dt_min.strftime('%Y%m%d') + '\n')
            fout.write('/end_date=' + dt_max.strftime('%Y%m%d') + '\n')
            fout.write('/start_time=' + dt_min.strftime('%H:%M:%S') + '[GMT]\n')
            fout.write('/end_time=' + dt_max.strftime('%H:%M:%S') + '[GMT]\n')

            #remove dt from out dict
            data_out.pop('dt', None)

            #output lat/lon headers
            lat_min = min(lat_lis)
            lat_max = max(lat_lis)
            fout.write('/north_latitude={:.3f}[DEG]\n'.format(lat_max))
            fout.write('/south_latitude={:.3f}[DEG]\n'.format(lat_min))

            lon_min = min(lon_lis)
            lon_max = max(lon_lis)
            fout.write('/east_longitude={:.3f}[DEG]\n'.format(lon_max))
            fout.write('/west_longitude={:.3f}[DEG]\n'.format(lon_min))

            #output remaining headers
            fout.write('/missing=' + missing + '\n')
            fout.write('/delimiter=comma\n')
            delim = ','

            #output comment headers
            dt_now = datetime.now()
            fout.write('! created ' + dt_now.strftime('%d %b %Y') + ' using awr2env.py\n')

            #output fields header
            fout.write('/fields=' + ','.join(data_out.keys()) + '\n')
            fout.write('/units=' + ','.join(unit_out.values()) + '\n')
            fout.write('/end_header\n')

            #output data matrix
            for i in range(len(data_out['year'])):
                row_ls = []

                for field in data_out:

                    if is_number(data_out[field][i]):
                        if isnan(float(data_out[field][i])):
                            row_ls.append(missing)
                        else:
                            row_ls.append(str(data_out[field][i]))
                    else:
                        if 'nan' in data_out[field][i].lower():
                            row_ls.append(missing)
                        else:
                            row_ls.append(str(data_out[field][i]))

                fout.write(delim.join(row_ls) + '\n')

    return

def is_number(s):

            """
            is_number determines if a given string is a number or not, does not handle complex numbers
            returns True for int, float, or long numbers, else False
            syntax: is_number(str)
            """

            try:
                float(s) # handles int, long, and float, but not complex
            except ValueError:
                return False
            return True

if __name__ == "__main__": main()
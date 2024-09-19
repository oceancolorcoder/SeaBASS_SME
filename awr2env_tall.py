
"""
Process SeaBASS Above Water Radiometry (AWR) to .env or .env.all
D. Aurin, NASA/GSFC 2024-04-05

This is meant for "tall" SeaBASS files on station with wavelength as a column
"""

def main(dict_args):
    import argparse
    from pathlib import Path
    from datetime import datetime
    from math import isnan
    import re
    import csv
    from copy import copy
    from collections import OrderedDict
    import numpy as np
    import pytz
    from SB_support import readSB#, is_number

    parser = argparse.ArgumentParser()
    fields_req = ['rrs']
    fields_opt = ['bincount', 'es', 'ed', 'lw'] 
    fields_dep = ['depth', 'pressure'] # No input depth fields. Depth is 0.
    missing = '-9999'

    metadata = dict_args['metadata']
    print(f"Flag file: {dict_args['flag_file']}") # SeaBASS only or validation
    print(f"Datatype: {metadata['dataType']} Instrument: {metadata['instrument']} Subinstrument: {metadata['subInstrument']} ")
    if dict_args['all'] ==  True:
        print('NOMAD files')
    else:
        print('VALIDATION files')


    timezone = pytz.utc
    filein_flag = Path(dict_args['flag_file'])
    if filein_flag.exists():
        with open(filein_flag, 'r') as csvfile:
            csvtable = csv.reader(csvfile)
            data = list(csvtable)
            data.pop(0)
        data_array = np.array(data, dtype=np.uint16)
        flagDatetime = [timezone.localize(datetime(*x)) for x in data_array]
        flag = data_array[:,6]
    else:
        parser.error('ERROR: invalid --flag_file specified; does ' + filein_flag.name + ' exist?')

    nSamples = len(dict_args['seabass_file'])
    # Loop over input SB files
    for fIndx, fnamein in enumerate(dict_args['seabass_file']):
        # filein_sb = Path(fnamein.name) ### __main__ parser gives a different list type here.
        filein_sb = Path(fnamein)
        print(f'Input SeaBASS file: {filein_sb}')

        # read and verify SeaBASS file and required fields
        if filein_sb.exists():
            ds = readSB(filename=str(filein_sb), \
                        mask_missing=True, \
                        mask_above_detection_limit=True, \
                        mask_below_detection_limit=True, \
                        no_warn=True)
            # nSamples+=ds.length
        else:
            parser.error('ERROR: invalid --seabass_file specified; does ' + filein_sb.name + ' exist?')

        #check for valid cruise name
        if not 'cruise' in ds.headers:
            parser.error('ERROR: cruise name not found in ' + filein_sb.name)

        #check for valid fields
        if not 'fields' in ds.headers:
            parser.error('ERROR: fields not found in ' + filein_sb.name)

        #check for valid dtimes
        ds.dtime = ds.fd_datetime() # SB_support
        ds.dtime = [timezone.localize(dt) for dt in ds.dtime]

        if not ds.dtime:
            parser.error('ERROR: date-time not parsable in ' + filein_sb.name)

        if 'measurement_depth' in ds.headers:
            header_depth = True
            depth = float(ds.headers['measurement_depth'])
        else:
            depth = 0.0

        match = re.search('[DEG]',ds.headers['north_latitude'])
        lat =float(match.string[0:match.start()-1])         
        match = re.search('[DEG]',ds.headers['east_longitude'])
        lon =float(match.string[0:match.start()-1])   

        firstGoodIndx = 0
        if fIndx == 0:
            # define output vars
            out_dir = Path('./') # Could change this to be the same folder

        #    Cruise.DataType_Instrument_OptionalSubInstrument.FirstInvestigator.OptionalSubcruise.env 
            if dict_args['all'] == True:
                fileout_sb = \
                    f"{ds.headers['experiment']}_{ds.headers['cruise']}_{ds.pi.split('_')[1]}_AOP_{metadata['subInstrument']}.env.all"
                    # f"{ds.headers['cruise'].lower()}.{metadata['dataType'].lower()}_{metadata['instrument'].lower()}_{metadata['subInstrument'].lower()}.{ds.pi.lower()}.env.all"
            else:
                fileout_sb = \
                    f"{ds.headers['experiment']}_{ds.headers['cruise']}_{ds.pi.split('_')[1]}_AOP_{metadata['subInstrument']}.env"
                    # f"{ds.headers['cruise'].lower()}.{metadata['dataType'].lower()}_{metadata['instrument'].lower()}_{metadata['subInstrument'].lower()}.{ds.pi.lower()}.env"
            
            # lat_lis  = []
            # lon_lis  = []

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

            data_out['depth'] = []
            unit_out['depth'] = 'm'

            data_out['cloud'] = [] # Only represents cloud reported in field notes
            unit_out['cloud'] = '%'

            data_out['associated_files'] = []
            unit_out['associated_files'] = 'none'

            data_out['associated_file_types'] = []
            unit_out['associated_file_types'] = 'none'

            data_out['flags'] = []
            unit_out['flags'] = 'none'
        
        # Screen for flagged data
        # if difference in flagDatetime and ds.dtime[i] is within 10 seconds...
        # For "tall" files, all ds.dtime should be the same.
        dateTimeDiff = [ds.dtime[0] - fDt for fDt in flagDatetime]
        absDTdiffsec = [abs(x.total_seconds()) for x in dateTimeDiff]
        if min(absDTdiffsec) < 10:
            index = absDTdiffsec.index(min(absDTdiffsec))
            # print(f'Match found {absDTdiffsec[index]} seconds from flag file')
        else:
            print(f'No matching time found in flag file: {ds.dtime[0]}')

        if flag[index] != 0:
            if firstGoodIndx == 0:
                firstGoodIndx = fIndx
            binFlag = 2**0 + 2**16 + 2**17 # 0(AOP) 9(OBPG software) 14(Es) 15(Rrs) 16(hyper) 17(above-water) None for Lsky/Li/Lw, etc.

            header_depth = False
            depth_field  = ''
            fields_fou = []

            for var in ds.variables.keys():    
                
                #check for required fields
                # Set binFlag depending on dataset(s) available
                ## This will be reset with every file, so they all have to be the same ##
                if var == 'rrs':
                    binFlag = binFlag + 2**15 # 15(Rrs)
                elif var ==  'es' or var == 'ed':
                    binFlag = binFlag + 2**14 # 14(Es)

                for field in fields_req:
                    
                    ma = re.search('^' + field + '$', var)

                    if ma:
                        # fields_fou.append(field)
                        fields_fou.append(var)

                    #check for error fields
                    for esuf in ds.err_suffixes: # ['_cv', '_sd', '_se', '_bincount']
                        efield = field + esuf
                        ma = re.search('^' + efield + '$', var)

                        if ma:
                            fields_fou.append(efield)

                #check for optional fields
                for field in fields_opt:
                    ma = re.search('^' + field + '$', var)

                    if ma:
                        fields_fou.append(field)

                    #check for error fields
                    for esuf in ds.err_suffixes: # ['_cv', '_sd', '_se', '_bincount']
                        efield = field + esuf
                        ma = re.search('^' + efield + '$', var)

                        if ma:
                            fields_fou.append(efield)

            if not fields_fou:
                parser.error('ERROR: AWR data not found in ' + filein_sb.name)

            # for i in range(ds.length):                        

            #verify row has valid AWR data
            flag_nodat = 0

            for field in fields_fou:
                for val in ds.data[field]:
                    if isnan(val):
                        flag_nodat = flag_nodat + 1

                if flag_nodat == len(ds.data[field]):
                    parser.error('ERROR: AWR all nans found in ' + filein_sb.name)
                
            # Append data
            i = 0 # Same across a file
            data_out['dt'].append(ds.dtime[i])

            data_out['year'].append(ds.dtime[i].strftime('%Y'))
            data_out['month'].append(ds.dtime[i].strftime('%m'))
            data_out['day'].append(ds.dtime[i].strftime('%d'))
            data_out['hour'].append(ds.dtime[i].strftime('%H'))
            data_out['minute'].append(ds.dtime[i].strftime('%M'))
            data_out['second'].append(ds.dtime[i].strftime('%S'))

            # lat_lis.append(ds.data['lat'][i])
            # lon_lis.append(ds.data['lon'][i])

            data_out['lat'].append('{:.4f}'.format(lat))
            data_out['lon'].append('{:.4f}'.format(lon))

            data_out['depth'].append(depth)

            if 'cloud' in ds.data.keys():
                if not isnan(ds.data['cloud'][i]):
                    data_out['cloud'].append(ds.data['cloud'][i])

                else:
                    data_out['cloud'].append('nan')
            else:
                    data_out['cloud'].append('nan')

            data_out['associated_files'].append(filein_sb.name)
            data_out['associated_file_types'].append('env')

            data_out['flags'].append(binFlag)

            #fill required data into output dict
            # Tack on the wavelength to the data_out key, values
            wavelength = ds.data['wavelength']
            # Test that subsequent files follow the first file
            if fIndx == firstGoodIndx:
                firstWL = wavelength
            if wavelength != firstWL:
                parser.error('ERROR: Wavelengths do not match first file in ' + filein_sb.name)

            for field in fields_fou:
                for iWv, wv in enumerate(wavelength):
                    if '_' in field:
                        [f1,f2] = field.split('_')
                        fieldXXX = f1 + str(wv) + '_' + f2

                    else:
                        fieldXXX = field + str(wv)

                    if not fieldXXX in data_out:
                        data_out[fieldXXX] = []
                        unit_out[fieldXXX] = ds.variables[field][-1]

                    if not isnan(ds.data[field][iWv]):
                        data_out[fieldXXX].append(ds.data[field][iWv])

                    else:
                        data_out[fieldXXX].append('nan')

    #sort data_out by dtime
    for field in data_out:
        if not 'dt' in field:
            data_out[field] = [x for (y,x) in sorted(zip(data_out['dt'],data_out[field]), key=lambda pair: pair[0])]

    # create and fill fileout_sb
    print(f"{len(data_out['dt'])} records retained of {nSamples}")
    print('Creating', fileout_sb)

    with open(out_dir / fileout_sb, 'w') as fout:

        #output headers
        fout.write('/begin_header\n')
        fout.write('/cruise=' + ds.headers['cruise'].lower() + '\n')
        fout.write('/affiliations=' + ds.headers['affiliations'].lower() + '\n')
        fout.write('/investigators=' + ds.headers['investigators'].lower() + '\n')
        fout.write('/experiment=' + ds.headers['experiment'].lower() + '\n')
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
        lat_all = [float(lati) for lati in data_out['lat']]
        lat_min = min(lat_all)
        lat_max = max(lat_all)
        fout.write('/north_latitude={:.3f}[DEG]\n'.format(lat_max))
        fout.write('/south_latitude={:.3f}[DEG]\n'.format(lat_min))

        lon_all = [float(loni) for loni in data_out['lon']]
        lon_min = min(lon_all)
        lon_max = max(lon_all)
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

    return fileout_sb

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

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
        This program creates .env or env.all file(s) for a given SeaBASS data file(s)
        containing Above Water Radiometry (AWR) and outputting a standard set of
        headers and fields common to the ENV file format.

        ENV.all files are interim files that are then culled into ENV files which are
        used to define and load validation matchup targets and data into the SeaBASS
        matchup work-flow.

        REQUIRED inputs:
            1) --seabass_file=    a valid SeaBASS file(s) with latitude, longitude, and date-time, and one or more of:
                                Rrs, Es.
            2) --flag_file=       a csv file with time, flag (0, 1, 2 for reject, seabass-only, validation)
            3) --all=             True/[False] (True to write .env.all file including non-validation spectra)

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
    parser.add_argument('--flag_file', nargs='+', type=argparse.FileType('r'), required=True, help='''\
        REQUIRED: input flag file
        Must be a csv file containing date-time and flag.
        ''')
    parser.add_argument('--all', nargs='+', type=argparse.FileType('r'), required=True, help='''\
        # REQUIRED: file type designation
        # True for .env.all, False for .env
        # ''')

    args=parser.parse_args()
    dict_args=vars(args)
    main(dict_args)
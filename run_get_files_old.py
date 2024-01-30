# Environment: conda activate sbanalyst
import sys
sys.path.append('/Users/daurin/GitRepos/sbanalyst')
import get_files

iFile = './Mannino_VIIRS_validation_2019_foster.lis'

# This construct would require changing get_files to capture args

get_files.main(iFile)

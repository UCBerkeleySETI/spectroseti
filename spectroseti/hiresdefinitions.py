__author__ = 'nate'
#
# Contains some definitions (directories, filenames)
# specific to the HIRES system
#
# location of reduced spectra
spectra_dir_reduced = '/mir3/iodfitsdb'
# location of raw spectra
spectra_dir_raw = '/mir3/raw/'
# location to save image output
output_png_dir = '/mir/ew/ntellis/laser/output/las_images/'
# Location of iodine order correction coeffs
iodine_dir = '/media/nate/DATA/Spectra/iodine/'
# name of database
database_name = 'test_database'
# yinds filename
yinds_file = '/data/yinds.npy'
# positive directory
pos_png_directory = '/media/nate/Ubuntu Storage/data/positive_raws/'
# Default wavelength scale for HIRES reduced FITS data
hires_wavs = '/data/hires_default_wavelengths.npy'
# Targets to ignore in image generation and laser search
ignore_targets = ['horizon lock', 'unknown', 'DOME FLAT', '',
                  '-FOUL WEATHER-', 'test','EXCEPTION - NO TARGNAME',
                  'v', 'LRIS (N wrap)', 'Moon','MOSFIRE Flats',
                  'NIRC STEP 1', 'cass module', 'dome flat',
                  'f/15 top end', 'fcass STEP 3', 'tertiary module',
                  'zenith lock', 'FST request pos','zenith  lock',
                  'SEG EXCHG', 'NIRC STEP 2', 'f/25 top end',
                  'Pluto Solar1', 'north wrap', 'thar', 'FST request pos',
                  'fcass STEP 1', 'fcass STEP 2', 'special park',
                  'HORIZON STOW', 'Seg Xchg', 'SEGMENT EXCHG2', 'GRB']
emission_wavs = [6561.1224365,  4860.0909424,  4339.3667603,  4100.7015991,
                 3969.0740967,  3888.0728149,  3834.4223022,  3796.9458008,
                 3769.6856689,  3749.2129517,  3733.4332275,  3721.0070801,
                 3711.0421753,  3702.9266357,  3696.2271118,  3690.6317139,
                 3685.9091187,  3681.8865967, 3736., 3819.6, 3859.8, 3933.,
                 4026.2, 4226., 4471.4, 5875.5, 6717.]

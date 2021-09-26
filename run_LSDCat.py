import os

# run the LSDCat process
def run_LSDCat(cube_name,band,SN_threshold):
    
    # make SN_threshold a string
    SN_threshold = str(SN_threshold)
    
    # change directory to cube's location
    os.chdir('C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/'+cube_name+'/SDSS_'+band+'_band')
    
    # run first function
    os.system('python ../../../../../lsdcat/lsd_cc_spatial.py --gaussian -pc 0.5 -S 1 -N 2 -t 1 -i SDSS_'+band+
    '_residual_image_'+cube_name+'.fits -o '+cube_name+'_'+band+'_spatial.fits')
    
    # run second
    os.system('python ../../../../../lsdcat/lsd_cc_spectral.py -i '+cube_name+'_'+band+'_spatial.fits --FWHM=300 -o '+cube_name+'_'+band+'_spatial_spectral.fits')
    
    # run third
    os.system('python ../../../../../lsdcat/tools/s2n-cube.py -i '+cube_name+'_'+band+'_spatial_spectral.fits -o '+cube_name+'_'+band+'_spatial_spectral_sn.fits')
    
    # run fourth
    os.system('python ../../../../../lsdcat/lsd_cat.py -i '+cube_name+'_'+band+'_spatial_spectral_sn.fits -t '+SN_threshold+' -c '+cube_name+'_'+band+'_LSDCat.cat')
    
    # run fifth
    os.system('python ../../../../../lsdcat/lsd_cat_measure.py -ta '+SN_threshold+' -f SDSS_'+band+'_residual_image_'+cube_name+'.fits -ic '+cube_name+
    '_'+band+'_LSDCat.fits --rmin 3 -ff '+cube_name+'_'+band+'_spatial_spectral.fits -ffsn '+cube_name+'_'+band+'_spatial_spectral_sn.fits --fhdu DATA')
import numpy as np
import matplotlib.pyplot as plt
from photutils.aperture import EllipticalAperture
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse
from photutils.isophote import build_ellipse_model
from photutils.isophote import EllipseSample, Isophote, IsophoteList
from astropy.io import fits
import time

def fit_isophote(image_name: str, band: str, procedure: str = None, vmin = -3, vmax = 3, init_ellipse:(float,float,float,float,float) = (27, 27, 15, 0.15,30*np.pi/180)):
    '''
    This function fits isophotes to an image that has been passed through an SDSS filter
    
    Inputs:
    
        
    Example:
        fit_isophote('MAGPI1501_subcube','i',procedure='photom')
        
    '''
    
    #TODO add a check to see if band and procedure strings are valid
    
    # time starting run
    t_start = time.time()
    
    # path to band image
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/SDSS_"+band+"_band_of_"+image_name+".npy"
    
    data = np.load(file_name)
    
    
    # This code needs to run regardless and doesn't take long so can run each time
    # following tutorials from https://photutils.readthedocs.io/en/stable/isophote.html and https://mpdaf.readthedocs.io/en/latest/cube.html
    
    # create initial guess ellipse
    geometry = EllipseGeometry(x0=init_ellipse[0], y0=init_ellipse[1], sma=init_ellipse[2], eps=init_ellipse[3], pa=init_ellipse[4])
    
    # calculating and plotting the initial ellipse guess
    aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,geometry.sma*(1 - geometry.eps),geometry.pa)
    plt.imshow(data,origin = 'lower',vmin=vmin,vmax=vmax)
    aper.plot(color='red')
    plt.title('Isophote initial guess')
    plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/Isophote_initial_guess.png")
    
    # calculating isophotes from that initial guesss
    ellipse = Ellipse(data, geometry=geometry)
    isolist = ellipse.fit_image(sclip=3.0, nclip=3)
    
    # can then print numerous isophotes for different semi-major axes
    model_image = build_ellipse_model(data.shape, isolist)
    residual = data - model_image
    
    plt.figure()
    plt.imshow(data, origin = 'lower',vmin=vmin,vmax=vmax)
    
    # will plot the ellipses for given sma values
    smas = np.linspace(10,50,10)
    
    # also plot the largest ellipse
    smas = np.append(smas,max(isolist.sma))
    
    print('Plotting isophotes with sma of: \n',smas)
    
    # actually plotting the ellipses
    for sma in smas:
        iso = isolist.get_closest(sma)
        x, y, = iso.sampled_coordinates()
        plt.plot(x,y, color='white')
    plt.title('Resulting isophote model image')
    plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/Modelled_isophotes.png")
    
    plt.figure()
    plt.imshow(model_image, origin = 'lower',vmin=vmin,vmax=vmax)
    plt.title('Modelled isophotes')
    plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/model_image.png")
    
    plt.figure()
    plt.imshow(residual, origin = 'lower',vmin=vmin,vmax=vmax)
    plt.title('The noise remove by the modelling process')
    plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/noise_from_process.png")
    
    # this code takes a while and may not need every run so give user option
    if procedure == 'photom':
        # from Step 7 at https://github.com/astropy/photutils-datasets/blob/main/notebooks/isophote/isophote_example4.ipynb
        
        with fits.open("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+".fits") as hdu:
            cube = hdu[1].data
        
        # choosing a few wavelength slices to model first
        # qfits indexes from 1, so slice number 2000 in pyfits is 1999 in python
        wavelength_indexes = np.linspace(499,2999, 20)
        
        # if want to print all wavelengths
        #wavelength_indexes = range(len(cube))
        
        # append this in testing as this is wavelength of known lense for 1501
        wavelength_indexes = np.append(wavelength_indexes,1860)
        
        for wavelength_index in wavelength_indexes:
            isolist_model_ = []
            # convert wavelength index to int
            wavelength_index = round(wavelength_index)
            
            # print every 100th run
            #if wavelength_index%100 == 0:
            #    print('Running slice',wavelength_index)
            
            # TODO didn't work if just had [:], maybe because g.sma = 0 for this?
            for iso in isolist[1:]:
                # get the EllipseGeometry of each modelled ellipse
                g = iso.sample.geometry
                
                # using the model image so can subtract source
                sample = EllipseSample(cube[wavelength_index], g.sma, geometry=g, integrmode='median', sclip=3.0, nclip=3)
                sample.update()
                
                iso_ = Isophote(sample, 0, True, 0)
                isolist_model_.append(iso_)
                
                # constructing the isophote list from the result
                isolist_model = IsophoteList(isolist_model_)
                
            # then plot the intensity profile
            # only use isolist.intens[1:] as removed first elem of isolist_model
            plt.figure()
            
            # if want relative intensity, use this code
            # errors = isolist_model.intens/isolist.intens[1:] * np.sqrt((isolist_model.int_err / isolist_model.intens)**2 + (isolist.int_err[1:] / isolist.intens[1:])**2)
            # plt.errorbar(isolist_model.sma**0.25, (isolist_model.intens / isolist.intens[1:]),yerr=errors, fmt='o', markersize=4)
            
            # if just want intensity use
            # errors = isolist_model.intens * np.sqrt((isolist_model.int_err / isolist_model.intens)**2)
            # plt.errorbar(isolist_model.sma**0.25, (isolist_model.intens),yerr=errors, fmt='o', markersize=4)
            
            
            # minus the model intensity from the actual
            errors = isolist_model.intens/isolist.intens[1:] * np.sqrt((isolist_model.int_err / isolist_model.intens)**2 + (isolist.int_err[1:] / isolist.intens[1:])**2)
            plt.errorbar(isolist_model.sma**(1/4), abs(isolist.intens[1:]-isolist_model.intens),yerr=errors, fmt='o-', markersize=4)
            
            plt.title("Intensity Profile of wavelength slice "+str(wavelength_index+1))
            plt.xlabel('sma**(1/4)')
            plt.ylabel('Intensity difference between model and data')
            
            # can use ylim to enforce limits on y-axis to allow for better comparison
            #plt.ylim([0.4,1.3])
            # add wavelength+1 so that wavelength of file is that as shown in qfitsview
            plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/Intensity_profile_wavelength_slice_"+str(wavelength_index+1)+".png")
    t_end = time.time()
    print('Run completed in,',t_end-t_start,'seconds!')
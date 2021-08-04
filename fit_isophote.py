import numpy as np
import matplotlib.pyplot as plt
from photutils.aperture import EllipticalAperture
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse
from photutils.isophote import build_ellipse_model
from photutils.isophote import EllipseSample, Isophote, IsophoteList
from astropy.io import fits

def fit_isophote(image_name: str, band: str, procedure: str = None, vmin = -3, vmax = 3, init_ellipse:(float,float,float,float,float) = (27, 27, 15, 0.15,30*np.pi/180)):
    '''
    This function fits isophotes to an image that has been passed through an SDSS filter
    
    Inputs:
    
        
    Example:
        fit_isophote('MAGPI1201','i',procedure='photom')
        
    '''
    
    #TODO add a check to see if band and procedure strings are valid
    
    # path to band image
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/SDSS_"+band+"_band_of_"+image_name+".npy"
    
    data = np.load(file_name)
    
    
    # This code needs to run regardless and doesn't take long so can run each time
    # following tutorials from https://photutils.readthedocs.io/en/stable/isophote.html and https://mpdaf.readthedocs.io/en/latest/cube.html
    
    # create initial guess ellipse
    # TODO MAKE IT SO CAN ENTER THESE
    geometry = EllipseGeometry(x0=init_ellipse[0], y0=init_ellipse[1], sma=init_ellipse[2], eps=init_ellipse[3], pa=init_ellipse[4])
    
    # calculating and plotting the initial ellipse guess
    aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,geometry.sma*(1 - geometry.eps),geometry.pa)
    plt.imshow(data,origin = 'lower',vmin=vmin,vmax=vmax)
    aper.plot(color='red')
    plt.title('Isophote initial guess')
    plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/Isophote_initial_guess.png")
    # calculating isophotes from that initial guesss
    ellipse = Ellipse(data, geometry=geometry)
    isolist = ellipse.fit_image()
    
    # can then print numerous isophotes for different semi-major axes
    model_image = build_ellipse_model(data.shape, isolist)
    residual = data - model_image
    
    plt.figure()
    plt.imshow(data, origin = 'lower',vmin=vmin,vmax=vmax)
    
    # will plot the ellipses for given sma values
    smas = np.linspace(10,50,5)
    
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
    
    # This code takes a while and may not be needed every time
    if procedure == 'photom':
        # from Step 7 at https://github.com/astropy/photutils-datasets/blob/main/notebooks/isophote/isophote_example4.ipynb
        
        with fits.open("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+".fits") as hdu:
            cube = hdu[1].data
        
        # choosing a few wavelength slices to model first
        # qfits indexes from 1, so slice number 2000 in pyfits is 1999 in python
        wavelength_indexes = [499,1999,2999]
        
        for wavelength_index in wavelength_indexes:
            isolist_narrow_ = []
            
            # TODO didn't work if just had [:], maybe because g.sma = 0 for this?
            for iso in isolist[1:]:
                # get the EllipseGeometry of each modelled ellipse
                g = iso.sample.geometry
                
                # using the model image (just one wavelength slice for testing) so can subtract source
                sample = EllipseSample(cube[wavelength_index], g.sma, geometry=g, integrmode='median', sclip=3.0, nclip=3)
                sample.update()
                
                iso_ = Isophote(sample, 0, True, 0)
                isolist_narrow_.append(iso_)
                
                # constructing the isophote list from the result
                isolist_narrow = IsophoteList(isolist_narrow_)
                
            # then plot the intensity profile
            # only use isolist.intens[1:] as removed first elem of isolist_narrow
            # TODO make it so don't have to do this??
            plt.figure()
            errors = isolist_narrow.intens/isolist.intens[1:] * np.sqrt((isolist_narrow.int_err / isolist_narrow.intens)**2 + (isolist.int_err[1:] / isolist.intens[1:])**2)
            plt.errorbar(isolist_narrow.sma**0.25, (isolist_narrow.intens / isolist.intens[1:]),yerr=errors, fmt='o', markersize=4)
            plt.title("Intensity Profile of wavelength slice "+str(wavelength_index+1))
            plt.xlabel('sma**1/4')
            plt.ylabel('Ratio')
            # add wavelength+1 so that wavelength of file is that as shown in qfitsview
            plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/"+image_name+"/SDSS_"+band+"_band/Intensity_profile_wavelength_slice_"+str(wavelength_index+1)+".png")
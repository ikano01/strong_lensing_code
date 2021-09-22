from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def filter_im(image_name: str, plot_single: bool = False, band: str = None, mask: bool = False):
    '''
    This function creates broadband images 
    and saves the data from each filter as .npy files
    
    Inputs:
        image_name: str
            A string containing the file name of the data cube fits file
        plot_single: bool
            Default is False, if True, plots a single image 
            in the band specified by the band variable
        band: str
            Default is None, but if want to plot a single image, must specify 
            what image band to plot. Can be 'g', 'r', 'i' or 'z'
        mask: bool
            Default is False, but if mask == True and a mask is defined in 
            the code that region will be masked
    
    Output:
        If plot_single == False:
            For each band, saves the intensity of each pixel as a numpy array
            Also saves a broadband image showing each filter
        If plot_single == True:
            Does same as above, but also saves a plot of a certain band
    
    Example:
        filter_im('MAGPI1201',plot_single = True, band = 'g',mask = True)
        Which saves numpy arrays of the intensity of each pixel for each band,
        saves a broadband image of all bands 
        and saves an image of just the 'g' band intensities
        
        If want to create image of larger subcube can do something like
        filter_im('MAGPI1501_larger_subcube')
    '''
    
    # loading in data cube
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
data_cubes/"+image_name+".fits"
    obj1 = Cube(file_name)
    
    # masking out any galaxies or areas specified
    if mask == True:
        # center = centre of mask in pixel units
        # radius = radius of mask in pixel units
        # if float, get a cicle mask but if (float,float) get rectangle mask

        obj1.mask_region(center= (15,36), radius = 2, unit_center = None, unit_radius = None)
    
    # getting images of different bands
    im_g = obj1.get_band_image('SDSS_g')
    im_r = obj1.get_band_image('SDSS_r')
    im_i = obj1.get_band_image('SDSS_i')
    im_z = obj1.get_band_image('SDSS_z')
    
    # saving the data out of each band as numpy arrays
    # each element is the intenstity of that pixel
    save_path = 'C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/'+image_name
    
    band_paths = [
        save_path+'/SDSS_g_band/',
        save_path+'/SDSS_r_band/',
        save_path+'/SDSS_i_band/',
        save_path+'/SDSS_z_band/'
        ]
    
    # creating the save directory if it doesn't exist
    for path in band_paths:
        if os.path.exists(path) == False:
            os.makedirs(path)
        
    # saving each file
    np.save(save_path+'/SDSS_g_band/SDSS_g_band_of_'+
            image_name+'.npy',ma.getdata(im_g))
    np.save(save_path+'/SDSS_r_band/SDSS_r_band_of_'+
            image_name+'.npy',ma.getdata(im_r))
    np.save(save_path+'/SDSS_i_band/SDSS_i_band_of_'+
            image_name+'.npy',ma.getdata(im_i))
    np.save(save_path+'/SDSS_z_band/SDSS_z_band_of_'+
            image_name+'.npy',ma.getdata(im_z))
    
    
    # plotting each
    fig, axs = plt.subplots(2, 2, figsize=(7, 7))
    
    im_g.plot(ax=axs[0,0], title='SDSS_g',zscale = True)
    im_r.plot(ax=axs[0,1], title='SDSS_r',zscale = True)
    im_i.plot(ax=axs[1,0], title='SDSS_i',zscale = True)
    im_z.plot(ax=axs[1,1], title='SDSS_z',zscale = True)
    
    plt.savefig(save_path+'/broadband_im.png')
    plt.show()
    
    if plot_single == True:
        if band not in ['g','r','i','z']:
            print('Cannot create single image, must specify band g,r,i or z')
            return
        
        # to plot single image
        # indexing the desired band image local variable to plot
        # eg. locals()['im_'+'g'] gives the variable im_g
        plot_im = locals()['im_'+band]
        plot_im.plot(title='SDSS_'+band,zscale = True)
        plt.savefig(save_path+'/SDSS_'+band+'_band/filter_image.png')
from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def filter_im(image_name: str, band: str = None, plot_single: bool = False):
    '''
    This function creates broadband images and then saves the data from each filter as .npy files
    
    Example:
        filter_im('MAGPI1201',band = 'g',plot_single = True)
    '''
    # loading in data cube
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+".fits"
    obj1 = Cube(file_name)
    
    # getting images of different bands
    im_g = obj1.get_band_image('SDSS_g')
    im_r = obj1.get_band_image('SDSS_r')
    im_i = obj1.get_band_image('SDSS_i')
    im_z = obj1.get_band_image('SDSS_z')
    
    # saving the data out of each band as numpy arrays
    save_path = 'C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/python files/data/'+image_name
    
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
    np.save(save_path+'/SDSS_g_band/SDSS_g_band_of_'+image_name+'.npy',ma.getdata(im_g))
    np.save(save_path+'/SDSS_r_band/SDSS_r_band_of_'+image_name+'.npy',ma.getdata(im_r))
    np.save(save_path+'/SDSS_i_band/SDSS_i_band_of_'+image_name+'.npy',ma.getdata(im_i))
    np.save(save_path+'/SDSS_z_band/SDSS_z_band_of_'+image_name+'.npy',ma.getdata(im_z))
    
    
    # plotting each
    fig, axs = plt.subplots(2, 2, figsize=(15, 7))
    
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
        # globals() converts what is in [] to be a variable
        # eg. globals()['im_'+'g'] == im_g
        plot_im = globals()['im_'+band]
        
        plot_im.plot(title='SDSS_'+band,zscale = True, cmap='gray')
        plt.savefig(save_path+'/SDSS_'+band+'_band/filter_image.png')
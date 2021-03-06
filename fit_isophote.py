import numpy as np
import matplotlib.pyplot as plt
from photutils.aperture import EllipticalAperture
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse
from photutils.isophote import build_ellipse_model
from photutils.isophote import EllipseSample, Isophote, IsophoteList
from astropy.io import fits
from mpdaf.obj import Cube
import time
import scipy

def fit_isophote(subcube_image_name: str, 
                 band: str, 
                 vmin = -3, vmax = 3, 
                 init_ellipse:[float,float,float,float,float]
                 = (27, 27, 15, 0.15,30), 
                 fit_isophote_guess: bool = False,
                 sclip_val=3.0, 
                 nclip_val=2, 
                 step_val=0.1,
                 larger_subcube = False):
    '''
    This function fits isophotes to an image that has been passed through an 
    SDSS filter.
    
    If larger_subcube == True, the isophote ellipses will be built using the larger subcube created earlier.
    All of the isophote images such as the initial isophote guess etc. will be saved in the larger subcube directory.
    
    The isophote ellipses created from the larger subcube are then used by the smaller subcube to build the isophote models.
    
    Inputs:
    init_ellipse:[float,float,float,float,float] where init_ellipse[0] is the y-coordinate of the centre of the ellipse,
        init_ellipse[1] is the x-coordinate of the centre of the ellipse, init_ellipse[2] is the size of the semi-major axis,
        init_ellipse[3] is the ratio between the semi-major and minor axes and init_ellipse[4] is the angle the ellipse is on in degrees
    
    DEFINITION:
        In this script, if larger_subcube == True, the variable names will be as expected, 
        however if larger_subcube == False, any variable names with larger_ in them will actually relate to the one subcube 
        that has been input into the code
    
        
    Example:
        fit_isophote('MAGPI1501_subcube','i',fit_isophote_guess = True)
    or to create the isophote models from a larger subcube
        fit_isophote('MAGPI1501_subcube','i',init_ellipse = (42, 42, 15, 0.15,30),fit_isophote_guess = True, larger_subcube = True)
    or
        fit_isophote('MAGPI1201_subcube','r',init_ellipse = (30, 30, 15, 0.15,30),fit_isophote_guess = False)
        
    '''
    
    #TODO add a check to see if inputs are valid
    
    # time starting run
    t_start = time.time()
    
# =============================================================================
#     Creating isophote curves from data cube
#
#     Following tutorials from 
#     https://photutils.readthedocs.io/en/stable/isophote.html and 
#     https://mpdaf.readthedocs.io/en/latest/cube.html
#
#     if larger_subcube == True, want to open the larger subcube to create the model
#     and then save it to the smaller subcube
#     therefore, change the file_name and cube data to be the larger subcube
# =============================================================================
    
    
    # path to band image
    # so larger_subcube_image_name will be the name of the larger subcube if larger = True, else it will be the same as subcube_image_name
    if larger_subcube:
        # assuming subcube_image_name is 9 numbers or letters then _subcube eg. MAGPI1501_subcube
        # and the larger subcube is 9 numbers or letters then _larger_subcube eg. MAGPI1501_larger_subcube
        larger_subcube_image_name = subcube_image_name[:9]+'_larger'+subcube_image_name[9:]
    else:
        larger_subcube_image_name = subcube_image_name
    
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_subcube_image_name+"/SDSS_"+band+"_band/SDSS_"+band\
+"_band_of_"+larger_subcube_image_name+".npy"

    # if larger_subcube == False, larger_subcube_name is equivalent to the subcube's name
    larger_subcube_file_path = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
data_cubes/"+larger_subcube_image_name+".fits"
    
    # loading in the pixel intensities
    pixel_intensities = np.load(file_name)
    
    
    # opening data cube
    with fits.open(larger_subcube_file_path) as hdu:
            # if larger_subcube == False, larger_subcube is equivalent to the subcube's data
            larger_subcube_data = hdu[1].data
            # change the nans in the data to be 0 otherwise fitting didn't work
            larger_subcube_data = np.nan_to_num(larger_subcube_data)
    
    # create initial guess ellipse
    geometry = EllipseGeometry(x0=init_ellipse[1], y0=init_ellipse[0], 
                               sma=init_ellipse[2], eps=init_ellipse[3], 
                               pa=np.deg2rad(init_ellipse[4]))
    
    # calculating and plotting the initial ellipse guess
    aper = EllipticalAperture((geometry.x0, geometry.y0), 
                              geometry.sma,geometry.sma*(1 - geometry.eps),
                              geometry.pa)
    
    plt.figure()
    plt.imshow(pixel_intensities,origin = 'lower',vmin=vmin,vmax=vmax)
    aper.plot(color='red')
    plt.title('Isophote initial guess')
    
    # saving figure
    plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_subcube_image_name+"/SDSS_"+band+"_band/\
Isophote_initial_guess.png")


    if not fit_isophote_guess:
        # is fit_isophote_guess == False, just show the plots and nothing else
        plt.show()
    else:
    
        # calculating isophotes from that initial guesss
        ellipse = Ellipse(pixel_intensities, geometry=geometry)
        isolist = ellipse.fit_image(sclip=sclip_val, nclip=nclip_val,step=step_val)
        
        # print numerous isophotes for different semi-major axes
        # this creates a 2d array of the isophote model
        model_image = build_ellipse_model(pixel_intensities.shape, isolist)
        
        residual = pixel_intensities - model_image
        
        plt.figure()
        plt.imshow(pixel_intensities, origin = 'lower',vmin=vmin,vmax=vmax)
        
        # will plot ellipses for given sma values
        smas = np.linspace(5,200,100)
        
        # also include largest ellipse
        smas = np.append(smas,max(isolist.sma))
        
        print('Plotting isophotes with sma of: \n',smas)
        
        # plotting the isophote ellipses
        for sma in smas:
            iso = isolist.get_closest(sma)
            x, y, = iso.sampled_coordinates()
            plt.plot(x,y, color='white')
        
        plt.title('Resulting isophote model image')
        plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_subcube_image_name+"/SDSS_"+band+"_band/Modelled_isophotes.png")
        
        plt.figure()
        plt.imshow(model_image, origin = 'lower',vmin=vmin,vmax=vmax)
        plt.title('Modelled isophotes')
        plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_subcube_image_name+"/SDSS_"+band+"_band/model_image.png")
        
        plt.figure()
        plt.imshow(residual, origin = 'lower',vmin=vmin,vmax=vmax)
        plt.title('The noise removed by the modelling process')

        plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_subcube_image_name+"/SDSS_"+band+"_band/noise_from_process.png")
        

    
# =============================================================================
#     create new model image from old isophote models but new intensities
# =============================================================================

        plt.figure()
        # just print a single slice
        plt.imshow(larger_subcube_data[2600], origin = 'lower',vmin=vmin,vmax=vmax)
        plt.title('A slice of the unblurred image')
        plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_subcube_image_name+"/SDSS_"+band+"_band/unblurred_image.png")

        # create copy of cube and apply a median filter for each pixel over 60 wavelength slices to further reduce noise
        larger_subcube_filtered = scipy.ndimage.median_filter(larger_subcube_data, size=(60, 1,1))
        
        # then subtracting this filtered image from the larger subcube
        larger_subcube_data -= larger_subcube_filtered
        
        plt.figure()
        # just print a single slice
        plt.imshow(larger_subcube_data[2600], origin = 'lower',vmin=vmin,vmax=vmax)
        plt.title('A slice of the image with the median blur subtracted')
        plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_subcube_image_name+"/SDSS_"+band+"_band/median_blurred_image.png")
        
        # opening subcube as Cube type
        cube_Cube_format = Cube("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/\
data/data_cubes/"+subcube_image_name+".fits")
        
        # create new empty cube to store residual wavelength slices
        # data_init and var_init initialise hdu[1] and hdu[2] to be filled
        residual_cube = cube_Cube_format.clone(data_init = np.zeros, var_init = np.zeros)
        
        # if want to print all wavelengths
        wavelength_indexes = range(len(larger_subcube_data))
        #wavelength_indexes = [1860]
        
        # if larger_subcube = True, also need the data from the smaller 
        #subcube so can subtract the model from as do not have it yet
        if larger_subcube:
            with fits.open("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/\
data/data_cubes/"+subcube_image_name+".fits") as hdu:
                subcube_data = hdu[1].data
                
                # need to also blur the smaller subcube
                subcube_filtered = scipy.ndimage.median_filter(subcube_data, size=(60, 1,1))
                subcube_data -= subcube_filtered
            
            # also need to find distance to crop the larger subcube
            # distance to crop larger cube to make it same as smaller cube
            # want length of just one cube slice
            dist_in = int(np.floor((len(larger_subcube_data[0])-len(subcube_data[0]))/2))
        
        for wavelength_index in wavelength_indexes:
            
            # print every 500th, 1000th etc. wavelength index
            if wavelength_index%100 == 0:
                print('Running wavelength slice '+str(wavelength_index))
            
            isolist_model_ = []
            # convert wavelength index to int
            wavelength_index = int(wavelength_index)
            
            # TODO didn't work if just had [:], maybe g.sma = 0 for this?
            for iso in isolist[1:]:
                # get the EllipseGeometry of each modelled ellipse
                g = iso.sample.geometry
                # using the model image so can subtract source
                sample = EllipseSample(larger_subcube_data[wavelength_index], g.sma, 
                                       geometry=g, integrmode='median', 
                                       sclip=sclip_val, nclip=nclip_val,astep=step_val)
                
                sample.update()
                
                iso_ = Isophote(sample, 0, True, 0)
                isolist_model_.append(iso_)
                
                # constructing the isophote list from the result
                isolist_model = IsophoteList(isolist_model_)
            
            # get the shape of the wavelength slice from the larger subcube if larger_subcube = True, else just load data of subcube
            new_slice_model_image = build_ellipse_model(larger_subcube_data[wavelength_index].shape,isolist_model)
            
            # creaing residual image for slice and saving it in blank cube
            if larger_subcube:
                # if the subcube_data slice and the subcube_data slice are of opposite parity, need to crop one side in one pixel more than the other
                if len(subcube_data[0])%2 != len(subcube_data[0])%2:
                    # [dist_in:(len(larger_subcube_data[0])-dist_in),dist_in:(len(larger_subcube_data[0])-dist_in)] is so crop model of larger subcube to same size as subcube
                    new_slice_model_image = new_slice_model_image[dist_in:(len(larger_subcube_data[0])-(dist_in+1)),dist_in:(len(larger_subcube_data[0])-(dist_in+1))]
                else:
                    # but if array even length, just crop using dist_in
                    new_slice_model_image = new_slice_model_image[dist_in:(len(larger_subcube_data[0])-dist_in),dist_in:(len(larger_subcube_data[0])-dist_in)]
            
                residual_cube[wavelength_index] = subcube_data[wavelength_index] - new_slice_model_image
            else:
                residual_cube[wavelength_index] = larger_subcube_data[wavelength_index] - new_slice_model_image
            
        # save variance array from original smaller subcube in residual subcube
        with fits.open("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
data_cubes/"+subcube_image_name+".fits") as hdu:
                
                # this is used when creating residual_cube
                variance = hdu[2].data
        
        residual_cube.var = variance
        
        # saving the cube
        new_save_path = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/\
data/python files/data/"+subcube_image_name+"/SDSS_"+band+"_band/SDSS_"+band+\
"_residual_image_"+subcube_image_name+".fits"
        residual_cube.write(new_save_path)
        print("Saved cube at",new_save_path)
        
        t_end = time.time()
        print('Run completed in,',t_end-t_start,'seconds!')
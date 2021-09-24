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

def fit_isophote(image_name: str, 
                 band: str, 
                 vmin = -3, vmax = 3, 
                 init_ellipse:(float,float,float,float,float) 
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
    # so larger_cube_image_name will be the name of the larger subcube if larger = True, else it will be the same as image_name
    if larger_subcube:
        # assuming image_name is 9 numbers or letters then _subcube eg. MAGPI1501_subcube
        # and the larger subcube is 9 numbers or letters then _larger_subcube eg. MAGPI1501_larger_subcube
        larger_cube_image_name = image_name[:9]+'_larger'+image_name[9:]
    else:
        larger_cube_image_name = image_name
    
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_cube_image_name+"/SDSS_"+band+"_band/SDSS_"+band\
+"_band_of_"+larger_cube_image_name+".npy"


    cube_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
data_cubes/"+larger_cube_image_name+".fits"
    
    # loading in the pixel intensities
    data = np.load(file_name)
    
    # opening data cube
    with fits.open(cube_name) as hdu:
            
            cube = hdu[1].data
    
    # create initial guess ellipse
    geometry = EllipseGeometry(x0=init_ellipse[0], y0=init_ellipse[1], 
                               sma=init_ellipse[2], eps=init_ellipse[3], 
                               pa=np.deg2rad(init_ellipse[4]))
    
    # calculating and plotting the initial ellipse guess
    aper = EllipticalAperture((geometry.x0, geometry.y0), 
                              geometry.sma,geometry.sma*(1 - geometry.eps),
                              geometry.pa)
    
    plt.figure()
    plt.imshow(data,origin = 'lower',vmin=vmin,vmax=vmax)
    aper.plot(color='red')
    plt.title('Isophote initial guess')
    
    # saving figure
    plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/"+larger_cube_image_name+"/SDSS_"+band+"_band/\
Isophote_initial_guess.png")


    if not fit_isophote_guess:
        # is fit_isophote_guess == False, just show the plots and nothing else
        plt.show()
    else:
    # if fit_isophote_guess == True, don't show plots but continue
    
        # calculating isophotes from that initial guesss
        ellipse = Ellipse(data, geometry=geometry)
        isolist = ellipse.fit_image(sclip=sclip_val, nclip=nclip_val,step=step_val)
        
        # print numerous isophotes for different semi-major axes
        # this creates a 2d array of the isophote model
        model_image = build_ellipse_model(data.shape, isolist)
        
        residual = data - model_image
        
        plt.figure()
        plt.imshow(data, origin = 'lower',vmin=vmin,vmax=vmax)
        
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
    python files/data/"+larger_cube_image_name+"/SDSS_"+band+"_band/Modelled_isophotes.png")
        
        plt.figure()
        plt.imshow(model_image, origin = 'lower',vmin=vmin,vmax=vmax)
        plt.title('Modelled isophotes')
        plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
    python files/data/"+larger_cube_image_name+"/SDSS_"+band+"_band/model_image.png")
        
        plt.figure()
        plt.imshow(residual, origin = 'lower',vmin=vmin,vmax=vmax)
        plt.title('The noise removed by the modelling process')

        plt.savefig("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
    python files/data/"+larger_cube_image_name+"/SDSS_"+band+"_band/noise_from_process.png")
        

    
# =============================================================================
#     create new model image from old isophote models but new intensities
# =============================================================================
        
        # if want to print all wavelengths
        wavelength_indexes = range(len(cube))
        
        # opening subcube as Cube type
        cube_Cube_format = Cube("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/\
data/data_cubes/"+image_name+".fits")
        
        # create new empty cube to store residual wavelength slices
        # data_init and var_init initialise hdu[1] and hdu[2] to be filled
        residual_cube = cube_Cube_format.clone(data_init = np.zeros, var_init = np.zeros)
        
        # if larger_subcube = True, also need the data from the smaller 
        #subcube so can subtract the mode from as do not have it yet
        if larger_subcube:
            with fits.open("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/\
data/data_cubes/"+image_name+".fits") as hdu:
                old_cube = hdu[1].data
        
        for wavelength_index in wavelength_indexes:
            isolist_model_ = []
            # convert wavelength index to int
            wavelength_index = round(wavelength_index)
            
            # if larger_subcube == True also need to find distance to crop the larger subcube
            if larger_subcube:
                # distance to crop larger cube to make it same as smaller cube
                # want length of just one cube slice
                dist_in = round((len(cube[0])-len(old_cube[0]))/2)
            
            # TODO didn't work if just had [:], maybe g.sma = 0 for this?
            for iso in isolist[1:]:
                # get the EllipseGeometry of each modelled ellipse
                g = iso.sample.geometry
                # using the model image so can subtract source
                sample = EllipseSample(cube[wavelength_index], g.sma, 
                                       geometry=g, integrmode='median', 
                                       sclip=sclip_val, nclip=nclip_val,astep=step_val)
                
                sample.update()
                
                iso_ = Isophote(sample, 0, True, 0)
                isolist_model_.append(iso_)
                
                # constructing the isophote list from the result
                isolist_model = IsophoteList(isolist_model_)
            
            # get the shape of the wavelength slice from the larger subcube if larger_subcube = True, else just load data of subcube
            new_slice_model_image = build_ellipse_model(cube[wavelength_index].shape,isolist_model)
            
            # creaing residual image for slice and saving it in blank cube
            if larger_subcube:
                # [dist_in:(len(cube[0])-dist_in),dist_in:(len(cube[0])-dist_in)] is so crop model of larger subcube to same size as subcube
                new_slice_model_image = new_slice_model_image[dist_in:(len(cube[0])-dist_in),dist_in:(len(cube[0])-dist_in)]
            
                residual_cube[wavelength_index] = old_cube[wavelength_index] - new_slice_model_image
            else:
                residual_cube[wavelength_index] = cube[wavelength_index] - new_slice_model_image
            
        # save variance array from original smaller subcube in residual subcube
        with fits.open("C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
data_cubes/"+image_name+".fits") as hdu:
                
                # this is used when creating residual_cube
                variance = hdu[2].data
        
        residual_cube.var = variance
        
        # saving the cube
        new_save_path = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/\
data/python files/data/"+image_name+"/SDSS_"+band+"_band/SDSS_"+band+\
"_residual_image_"+image_name+".fits"
        residual_cube.write(new_save_path)
        print("Saved cube at",new_save_path)
        
        t_end = time.time()
        print('Run completed in,',t_end-t_start,'seconds!')
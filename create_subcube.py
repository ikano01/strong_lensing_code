from mpdaf.obj import Cube
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def create_subcube(image_name:str, centre:(int,int),radius:int):
    '''
    This function will create and save a new subcube.
    
    Example:
        create_subcube('MAGPI1201', (197,197), 25)
    '''
    # TODO Need to fix this file!!! Does not work


    # load in full cube
    file_name = "C:/Users/isaac/Documents/Uni 2021/Sem 2/ASTR3005/data/data_cubes/"+image_name+".fits"
    
    with fits.open(file_name) as hdul:
        hdul.info()
        
        # opening DATA
        data = hdul[1].data
        print('\n Loaded data of shape',np.shape(data))
    
    # define coordinates to cut along
    x_coord_subset = [centre[0]-radius,centre[0]+radius]
    y_coord_subset = [centre[1]-radius,centre[1]+radius]
    
    # create subcube
    data = np.array(data)
    subcube = data[:,y_coord_subset,x_coord_subset]
    
    # saving subcube
    save_name = "C:/Users/isaac/Documents/Uni 2021/Sem 2/ASTR3005/data/data_cubes/"+image_name+"_subcube.fits"


create_subcube('MAGPI1201', (197,197), 25)
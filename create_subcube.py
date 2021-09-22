from mpdaf.obj import Cube
from mpdaf.obj import WCS
import astropy.units as u
import matplotlib.pyplot as plt

def create_subcube(image_name:str, centre:[float,float], length_coord:[float,float], subtype:str = 'cube', length_axis:str = 'x', larger_subcube:str = False):
    '''
    This function will create and save a new subcube.
    The option for larger_subcube = True can be set, which just changes the name of the output file so as to not affect previous smaller subcubes.
    
    MPDAF images are stored in python arrays that are indexed in [y,x] axis order
    
    Examples:
        create_subcube('MAGPI1201', [197,197], [256,197],subtype = 'cube',length_axis='y')
    or
        create_subcube('MAGPI1501', [210,192], [263,192],subtype = 'cube',length_axis='y')
    or to create larger subcube cut
        create_subcube('MAGPI1501', [210,192], [293,192],subtype = 'cube',length_axis='y', larger_subcube = True)
    or
        create_subcube('MAGPI1203', [197,197], [260,197],subtype = 'cube',length_axis='y')
    or
        create_subcube('MAGPI1508', [195,200], [235,195],subtype = 'cube',length_axis='y')
    '''
    
    # load in full cube
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+".fits"
    full_cube = Cube(file_name)    
    
    
    ### creating subcube ###
    
    # Due qfitsview to indexing from 1, need to subtract 1 from the pixel count
    centre[0] -= 1
    centre[1] -= 1
    length_coord[0] -= 1
    length_coord[1] -= 1
    
    # Instancing the WCS class, which deals with world coordinates
    wcs = WCS(full_cube.get_wcs_header())
    # converting centre into degrees
    centre = wcs.pix2sky(centre)[0]
    print('Centre in (ra,dec):',centre)
    
    
    # converting length coordinate to degrees
    # TODO MIGHT NEED TO INPUT BOX LENGTH AS COORDINATE SO HAVE Y AND X AND THEN ONLY INPUT ONE OF THEM INTO THE SUBCUBE FUNCTION
    # had to index [0] because wcs.pix2sky output the coordinates in a second list, eg. [[y,x]]
    length_coord = wcs.pix2sky(length_coord)[0]
    
    # getting the length in degrees by subtracting y-degree (or x) of centre with y-degree (or x) of length_coord
    if length_axis == 'x':
        length = abs(centre[1]-length_coord[1])
    elif length_axis == 'y':
        length = abs(centre[0]-length_coord[0])
    else:
        raise SyntaxError('Must specify whether to use the y or x value of the length coordinate [y,x]')
    print('Cube side length in degrees:',length)
    
    # converting length to arcseconds
    # .value just gives the number without the unit
    length = ((length*u.deg).to(u.arcsec)).value
    print('Cube side length in arcseconds:',length)
    
    
    if subtype == 'cube':
        # creating cubic subcube
        subcube = full_cube.subcube(center = (centre[0],centre[1]),size = length)
    elif subtype == 'circle':
        # creating circle aperature subcube
        subcube = full_cube.subcube_circle_aperture(centre,radius=length)
    else:
        raise SyntaxError('Must specify sub_type of either cube or circle')
        return    
    
    # saving subcube
    if larger_subcube:
        save_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+"_larger_subcube.fits"
    else:
        save_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+"_subcube.fits"
    
    subcube.sum(axis=0).plot()
    plt.show()
    subcube.write(save_name)


create_subcube('MAGPI1501', [210,192], [293,192],subtype = 'cube',length_axis='y')
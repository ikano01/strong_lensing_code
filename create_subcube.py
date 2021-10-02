from mpdaf.obj import Cube
from mpdaf.obj import WCS
import astropy.units as u
import matplotlib.pyplot as plt

def create_subcube(image_name:str, centre:[int,int], length_coord:[int,int], subtype:str = 'cube', length_axis:str = 'x', 
    mask_cube:bool = False, larger_subcube:bool = False, pixels_larger:int = 0, mask_centres:list = [], mask_extents:list = []):
    '''
    This function will create and save a new subcube.
    
    Inputs:
        cube_name:str is the name of the downloaded full cube
        centre:[int,int] is the coordinate in terms of [y,x] that will be the centre of the subcube
        length_coord:[int,int] is the coordinate of the distance to cut to, either vertically or horizontally from the centre_coord depending if length_axis is 'x' or 'y'
        length_axis:str can be either 'y' if the length_coord is vertical from the centre_coord or 'x' if the length_coord is horizontal from the centre_coord
        mask_cube:bool if True adds a mask to the subcube
        mask_centres:list is a list of lists which each correspond to the [y,x] centre of a mask
        mask_extents:list is a list of lists which each correspond to the [y,x] extents of the mask
        larger_subcube:bool if True, this will save the subcube with _larger added to the filename
        pixels_larger:int is the number of pixels wider and higher the new subcube should be
    
    MPDAF images are stored in python arrays that are indexed in [y,x] axis order
    
    Examples:
        create_subcube('MAGPI1201', [197,197], [256,197],subtype = 'cube',length_axis='y')
    or
        create_subcube('MAGPI1501', [210,192], [263,192],subtype = 'cube',length_axis='y')
    or to create larger subcube cut
        create_subcube('MAGPI1501', [210,192], [293,192],subtype = 'cube',length_axis='y', larger_subcube = True, pixels_larger = 30)
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
        if larger_subcube:
            # if larger_subcube, want to make length pixels_larger bigger
            length = abs(centre[1]-length_coord[1])+pixels_larger
        else:
            length = abs(centre[1]-length_coord[1])
    elif length_axis == 'y':
        # if larger_subcube, want to make length pixels_larger bigger
        if larger_subcube:
            length = abs(centre[0]-length_coord[0])+pixels_larger
        else:
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
    
    # masking the cube
    ### TODO won't actually be able to see the subcube when set these so need to be able to put in centre and extents etc. to be for the full subcube
    ### maybe could translate it by number of pixels away from the subcube centre???
    if mask_cube:
        for i,mask_centre in enumerate(mask_centres):
            # mask extents actually have to be input as (x,y), though I thought this would be confusing to input as the rest is (y,x) so need to change it here
            mask_extent = [mask_extents[i][1],mask_extents[i][0]]
            
            if larger_subcube:
                # to fit the mask to the correct position in larger subcube, needs to be shifted by pixels_larger/2
                mask_centre = [mask_centre[0]+pixels_larger/2,mask_centre[1]+pixels_larger/2]
            
            # I set the rotation to 0, but if need it can implement it by lettting user change 0 to something else
            # unit_center and unit_radius are set to None so the centre and extent values have units of pixels
            subcube.mask_ellipse(mask_centre,mask_extents,0,unit_center = None, unit_radius = None)
        pass
    
    
    # saving subcube
    if larger_subcube:
        save_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+"_larger_subcube.fits"
    else:
        save_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+"_subcube.fits"
    
    subcube.sum(axis=0).plot()
    plt.show()
    subcube.write(save_name,savemask='nan')
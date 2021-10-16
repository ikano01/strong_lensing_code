from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import gc

def create_subcube(image_name:str, centre:[int,int], length:[int,int], 
    mask_cube:bool = False, larger_subcube:bool = False, pixels_larger:int = 0, mask_centres:list = [], mask_extents:list = []):
    '''
    This function will create and save a new subcube.
    
    Inputs:
        cube_name:str is the name of the downloaded full cube
        centre:[int,int] is the coordinate in terms of [y,x] that will be the centre of the subcube
        length:int is the radius in pixels to cut to (so the whole cube has side lengths of 2*length)
        mask_cube:bool if True adds a mask to the subcube
        mask_centres:list is a list of lists which each correspond to the [y,x] centre of a mask to be placed on the full cube
        mask_extents:list is a list of lists which each correspond to the [y,x] radii of the mask to be placed on the full cube
        larger_subcube:bool if True, this will save the subcube with _larger added to the filename
        pixels_larger:int is the number of pixels wider and higher on each side the new subcube should be
    
    MPDAF images are stored in python arrays that are indexed in [y,x] axis order
    
    Examples:
        create_subcube('MAGPI1201', [197,197], 59)
    or
        create_subcube('MAGPI1501', [210,192], 53)
    or to create larger subcube cut
        create_subcube('MAGPI1501', [210,192], 83, larger_subcube = True, pixels_larger = 30)
    or
        create_subcube('MAGPI1203', [197,197], 63)
    or a mask example
        create_subcube('MAGPI1203', [197,197], 93,mask_cube = True, mask_centres = [[161,198]],mask_extents = [[9,9]])
    or
        create_subcube('MAGPI1508', [195,200], 40)
    '''
    
    # load in full cube
    file_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+".fits"
    full_cube = Cube(file_name)
    
    
    ### creating subcube ###
    
    # initialising centre_coord and length arrays
    centre_coord = [0,0]
    length_extent = [0,0]
    
    # Due qfitsview to indexing from 1, need to subtract 1 from the pixel count
    centre_coord[0] = centre[0] - 1
    centre_coord[1] = centre[1] - 1
    
    
    length_extent = length
    
    # if larger_subcube, want to make length pixels_larger bigger
    if larger_subcube:
        # want to make length_extent[1] pixels_larger further from the centre
        length_extent += pixels_larger
    
    # creating cubic subcube
    subcube = full_cube.subcube(center = centre_coord,size = 2*length_extent, unit_center = None, unit_size = None)
    
    # as full_cube no longer needed, delete to save space
    del full_cube
    # gc.collect() clears the variable from memory
    gc.collect()
    
    # masking the cube
    if mask_cube:
        for i,mask_centre in enumerate(mask_centres):
            # first need to convert the mask_centre coordinates from the full cube to the subcube

            # all of the pixels in the subcube shifted coordinates when cropped from the subcube.
            # as the subcube is a square centered on the centre pixel and centre pixel in full cube coordinates is known
            # the pixel shift is the difference between the centre full cube coord and its subcube coordinate which is 
            # at the centre of the cube, (so this is len(subcube.data[0])/2)
            # centre is the centre coords in the full cube coordinates
            # Also have to minus 1 because if len(subcube.data[0]) is 30, so len(subcube.data[0])/2 is 15, 
            # but the point actually has an index of 14
            ### Figure out why need to plus 1 - something to do with indexing from 1 ###
            if len(subcube.data[0])%2 == 0:
                pixel_shift = [centre_coord[0] - (len(subcube.data[0])/2-1),centre_coord[1] - (len(subcube.data[0])/2-1)]
            else:
                # if the width of the subcube is odd, must add one so that pixel shift is an integer
                # Also have to minus 1 (see above)
                pixel_shift = [centre_coord[0] - ((len(subcube.data[0])+1)/2-1),centre_coord[1] - ((len(subcube.data[0])+1)/2-1)]
            
            # initialising array
            shifted_mask_centre = [0,0]
            
            # the mask centre in the subcube is then the centre coordinate in the full cube coordinates minus the pixel_shift
            # must subtract by 1 because QFitsview indexes from 1
            shifted_mask_centre[0] = mask_centre[0] - pixel_shift[0] - 1
            shifted_mask_centre[1] = mask_centre[1] - pixel_shift[1] - 1
            
            # mask extents actually have to be input as (x,y), though I thought this would be confusing to input as the rest is (y,x) so need to change it here
            # also have to add 0.5 to each so the mask extends mask_extent pixels from the mask_centre
            mask_extent = [mask_extents[i][1]+0.5,mask_extents[i][0]+0.5]
            
            # I set the rotation to 0, but if need it can implement it by lettting user change 0 to something else
            # unit_center and unit_radius are set to None so the centre and extent values have units of pixels
            subcube.mask_ellipse(shifted_mask_centre,mask_extent,0,unit_center = None, unit_radius = None)
            print('Creating mask centred at [y,x]: ',shifted_mask_centre,' (when indexing from 0)')
    
    # saving subcube
    if larger_subcube:
        save_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+"_larger_subcube.fits"
    else:
        save_name = "C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/data_cubes/"+image_name+"_subcube.fits"
    
    # uncomment if want to see an example of the subcube
    #plt.figure()
    #subcube.sum(axis=0).plot()
    #plt.show()
    
    subcube.write(save_name,savemask='nan',convert_float32 = False)
    # delete the subcube to save space
    del subcube
    # gc.collect() clears the variable from memory
    gc.collect()
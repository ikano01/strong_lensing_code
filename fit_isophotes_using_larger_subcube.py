from create_subcube import create_subcube
from create_filters_image import filter_im
from fit_isophote import fit_isophote
from run_LSDCat import run_LSDCat
import numpy as np

### code to fit isophote to larger subcube and then copy onto a smaller cube ###

def run(cube_name,centre,length, mask_cube = False, mask_centres = [], mask_extents = [],
    pixels_larger = 30,init_ellipse = [30, 30, 15, 0.15,30],check_init_ellipse = True,band = 'g',SN_threshold = 10):
    '''
    Inputs:
        cube_name:str is the name of the downloaded full cube
        centre:[int,int] is the coordinate of the centre of the lens galaxy where coordinates are expressed [y,x]
        length:int is the radius in pixels to cut to (so the whole cube has side lengths of 2*length)
        mask_cube:bool if True adds a mask to the subcube
        mask_centres:list is a list of lists which each correspond to the [y,x] centre of a mask to be placed on the full cube
        mask_extents:list is a list of lists which each correspond to the [y,x] extents of the mask to be placed on the full cube - Can measure number of pixels using Qfitsview by pressing p
        pixels_larger:int is the number of pixels larger the larger subcube should be from the length_coord
        init_ellipse: (int,int,int,float,float) is the equestion for the initial ellipse guess for the larger subcube
        check_init_ellipse: bool if check_init_ellipse == True, the isophotes will not be subtracted straight away and instead 
            the initial ellipse given will be plotted on the larger subcube image so that the fit can be checked. The user can then decide to re-define the ellipse or 
            to continue running the simulation
        band:str is the band to run the isophote fit on
        SN_threshold: int is the signal to noise threshold that counts as a detection within the LSDCat processing steps
        
        eg.
        run('MAGPI1201',[197,197], 59, 30,[30, 30, 15, 0.15,30],check_init_ellipse = True,band='g', SN_threshold = 10)
        or
        run('MAGPI1203',[196,196], 65, 20,[45, 45, 15, 0.1,0],check_init_ellipse = True,band='g', SN_threshold = 10)
        or for the masked version,
        run(cube_name='MAGPI1203',centre=[196,196], length=65,mask_cube=True,mask_centres=[[153,197],[194,213],[238,191],[194,163],[153,195],[227,219],[205,245],[170,233],[
            170,249],[176,262]],mask_extents=[[24,35],[5,5],[17,25],[4,3],[11,7],[6,16],[10,7],[8,7],[7,7],[7,9]], pixels_larger=20,init_ellipse=[45, 45, 15, 0.1,0], 
            check_init_ellipse = True,band='g', SN_threshold = 10)
        or
        run(cube_name='MAGPI1205',centre=[191,199], length=38,mask_cube=True,mask_centres=[[157,197],[211,226],[223,177],[238,200]],
            mask_extents=[[13,16],[4,4],[7,7],[5,5]], pixels_larger=20,init_ellipse=[25, 25, 20, 0.1,0], 
            check_init_ellipse = True,band='g', SN_threshold = 10)
    '''
    # if the cube should be masked maskmask_cube == True
    # create subcube
    create_subcube(cube_name, centre, length,subtype = 'cube', larger_subcube = False, pixels_larger = pixels_larger, mask_cube = mask_cube, mask_centres = mask_centres,mask_extents = mask_extents)
    
    # create larger subcube
    create_subcube(cube_name, centre, length,subtype = 'cube', larger_subcube = True, pixels_larger = pixels_larger, mask_cube = mask_cube, mask_centres = mask_centres,mask_extents = mask_extents)

    # create filters image of both cubes
    filter_im(cube_name+'_subcube')
    filter_im(cube_name+'_larger_subcube')

    if check_init_ellipse:
        run_answer = 'N'
        while run_answer.upper() == 'N':
            print('The current initial ellipse model is:')
            fit_isophote(cube_name+'_subcube',band,init_ellipse = init_ellipse,fit_isophote_guess = False, larger_subcube = True)
            run_answer = input('Would you like to run the code with this ellipse? y or n? \n')
            if run_answer.upper() == 'N':
                # new_ellipse_form is then a list of the inputs
                print('The old form of the ellipse was: \n'+str(init_ellipse))
                for i,elem in enumerate(init_ellipse):
                    init_ellipse[i] = float(input('What would you like element '+str(i)+' to be?\nIt was originally '+str(init_ellipse[i])+'\n'))
        print('Running isophote fit!')
    
    # when the user is happy with the form, the code will then run the full fitting
    # fit the isophotes using the larger subcube as a base to build the models
    fit_isophote(cube_name+'_subcube',band,init_ellipse = init_ellipse,fit_isophote_guess = True, larger_subcube = True)
    
    # process the cube using the LSDCat files
    run_LSDCat(cube_name+'_subcube',band,SN_threshold)
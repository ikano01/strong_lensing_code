from create_subcube import create_subcube
from create_filters_image import filter_im
from fit_isophote import fit_isophote
import numpy as np

### code to fit isophote to larger subcube and then copy onto a smaller cube ###

def run(cube_name,centre_coord,length_coord,length_axis, pixels_larger,init_ellipse,check_init_ellipse,band):
    '''
    Inputs:
        cube_name:str is the name of the downloaded full cube
        centre_coord:[float,float] is the coordinate of the centre of the lens galaxy where coordinates are expressed [y,x]
        length_coord:[float,float] is the coordinate of the distance to cut to, either vertically or horizontally from the centre_coord depending if length_axis is 'x' or 'y'
        length_axis:str can be either 'y' if the length_coord is vertical from the centre_coord or 'x' if the length_coord is horizontal from the centre_coord
        pixels_larger:int is the number of pixels larger the larger subcube should be from the length_coord
        init_ellipse: (int,int,int,float,float) is the equestion for the initial ellipse guess for the larger subcube
        check_init_ellipse: bool if check_init_ellipse == True, the isophotes will not be subtracted straight away and instead 
            the initial ellipse given will be plotted on the larger subcube image so that the fit can be checked. The user can then decide to re-define the ellipse or 
            to continue running the simulation
        band:str is the band to run the isophote fit on
        
        eg.
        run('MAGPI1201',[197,197], [256,197],length_axis='y', [286,197],(30, 30, 15, 0.15,30*np.pi/180),check_init_ellipse = True)
    '''

    # create subcube
    create_subcube(cube_name, centre_coord, length_coord,subtype = 'cube',length_axis=length_axis)

    # create larger subcube
    if length_axis == 'x':
        larger_subcube_length_coord = [length_coord[0],length_coord[1]+pixels_larger[0]]
    else:
        larger_subcube_length_coord = [length_coord[1]+pixels_larger[1],length_coord[0]]
    
    create_subcube(cube_name, centre_coord, larger_subcube_length_coord,subtype = 'cube',length_axis=length_axis, larger_subcube = True)

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
    
    # TODO add in ability to run LSDCat commands eg.
    # os.system('python ../../../../../lsdcat/lsd_cc_spatial.py --gaussian -pc 0.5 -S 1 -N 2 -t 1 -i SDSS_i_residual_image_MAGPI1501_subcube.fits -o MAGPI1501_i_subcube_spatial.fits')

run('MAGPI1201',[197,197], [256,197],'y', [286,197],[30, 30, 15, 0.15,30],check_init_ellipse = True,band='r')
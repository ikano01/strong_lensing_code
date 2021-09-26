import ctypes

# load in the .dat file
file_name = 'C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/MAGPI1501_subcube/SDSS_i_band/\
SDSS_i_residual_image_MAGPI1501_subcube_MAGPI1501_i_subcube_LSDCat_IKA_class.dat'
file_name = 'C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/MAGPI1201_subcube/SDSS_r_band/\
SDSS_r_residual_image_MAGPI1201_subcube_MAGPI1201_r_subcube_LSDCat_IKA_class.dat'
file_name = 'C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/MAGPI1508_subcube/SDSS_r_band/\
SDSS_r_residual_image_MAGPI1508_subcube_MAGPI1508_r_subcube_LSDCat_IKA_class.dat'
file_name = 'C:/Users/isaac/Documents/Uni_2021/Sem_2/ASTR3005/data/\
python files/data/MAGPI1202_subcube/SDSS_r_band/\
SDSS_r_residual_image_MAGPI1202_subcube_MAGPI1202_subcube_r_LSDCat_IKA_class.dat'

### ONLY USE WITH COPY OF FILE JUST IN CASE ###

# reading in file
with open(file_name,'r') as file:
    file_lines = file.readlines()
    
for line_num,possible_source in enumerate(file_lines):
    # if line doesnt start as a comment - this gets rid of headers and comment lines
    # final number in line is at [-2]
    # possible_source[-2] means not galaxy
    if possible_source[0] != '#' and possible_source[-2] != '2':
        # assuming max pixel num is < 100
        # x elem is at [24:26]
        # y elem is at [33:35]
        # z elem is at [40:45]
        
        answer = ctypes.windll.user32.MessageBoxW(0, 'x: '+str(possible_source[24:26])+' y: '+str(possible_source[33:35])+' z: '+str(possible_source[40:44])+" Is this a galaxy \n Yes to make type 1, No for 2 and Cancel for 3", "Type "+str(possible_source[-2])+" Galaxy", 3)
        # on final line of dat file, it doesn't have a \n so need to do differently
        if line_num != len(file_lines)-1:
            # if click yes (returns 6) change classification to yes a galaxy
            if answer == 6:
                print('Yes')
                # go up to second last element as \n is counted as 1 character
                possible_source = possible_source[0:-2] + '1\n'
            
            # if click no (returns 7), change to not a galaxy
            elif answer == 7:
                print('No')
                possible_source = possible_source[0:-2] + '2\n'
            
            # if click cancel (returns 2), change to not sure
            elif answer == 2:
                print('Maybe')
                possible_source = possible_source[0:-2] + '3\n'
        else:
            # if click yes (returns 6) change classification to yes a galaxy
            if answer == 6:
                print('Yes')
                possible_source = possible_source[0:-1] + '1'
                # TODO now need to write back into file??
            
            # if click no (returns 7), change to not a galaxy
            elif answer == 7:
                print('No')
                possible_source = possible_source[0:-1] + '2'
            
            # if click cancel (returns 2), change to not sure
            elif answer == 2:
                print('Maybe')
                possible_source = possible_source[0:-1] + '3'
            
        
        
        file_lines[line_num] = possible_source
        print(possible_source)

print(file_lines)
# rewriting file
with open(file_name,'w') as file:
    file.writelines(file_lines)

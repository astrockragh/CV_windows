import numpy
import pickle
import scipy.optimize
import scipy.interpolate
import scipy.integrate


def linear(x,m,c):
    return m*x+c
#--------------------------------------------------------------------------------------------------------------
# This function takes the input parameters and computes the comving dimensions of the survey in units of Mpc/h
def get_survey_dimensions(central_redshift,survey_length_in_arc_sec,survey_aspect_ratio): # Calculates the survey dimensions in cMpc/h for a given field of view (in arc-sec) and redshift width 
    distance_to_central_redshift=comoving_distance(0,central_redshift)
    #print "comoving_distance:",distance_to_central_redshift/h
    survey_length_in_radians=survey_length_in_arc_sec*arc_sec_to_radians
    survey_length_in_cmpc_h=survey_length_in_radians*distance_to_central_redshift
    survey_dimensions=numpy.array([survey_length_in_cmpc_h,survey_length_in_cmpc_h*survey_aspect_ratio])
    return survey_dimensions
#---------------------------------------------------------------------------------------------------------------
# This function computes the comoving distance between two redshifts
def comoving_distance(z_low,z_up): # Calculates the comoving distance between two redshifts 
    box_size=scipy.integrate.quad(lambda z:c/H(z),z_low,z_up)[0]
    return box_size
#---------------------------------------------------------------------------------------------------------------

arc_sec_to_radians=1./(180/3.14)*(1./3600.)

f=open('./COSMIC_VARIANCE_AT_HIGHZ/cosmology.txt')
lines=f.readlines()
print("Cosmological parameters (Default: BlueTides cosmology. Edit 'cosmology.txt' to change)\n")
for line in lines:
    if ('-------' not in line):
        #line= line.replace('=',' ')
        #line=line.split()
        print(line)
        exec(line)
        
H = lambda z: H0*(om0*(1+z)**(3) + oml)**0.5


def get_cosmic_variance(H_band_threshold,redshift,redshift_width,survey_length_in_arc_sec,survey_aspect_ratio):
    import numpy
    import pickle
    import scipy.optimize
    import scipy.interpolate
    import scipy.integrate
    dimensions=get_survey_dimensions(redshift,survey_length_in_arc_sec,survey_aspect_ratio)
    length_in_Mpc=numpy.sqrt(dimensions[0]*dimensions[1]) 
    if (length_in_Mpc>400):
        print("Warning!!! Comoving dimensions of the survey is ",length_in_Mpc,"Mpc/h, greater than the BlueTides box length 400 Mpc/h. Do not trust the estimate")
    if (length_in_Mpc<0.01):
        print("Warning!!! Comoving dimensions of the survey is ",length_in_Mpc,"Mpc/h, less than the typical galaxy size 0.01 Mpc/h. Do not trust the estimate")
    
    reference_redshift_width=0.1
    redshift_space=[7.5,8.0,9.0,10.0]
    mUV_space=[25.,26.,28.0,29,30.0]
    sigma_space_redshift=[]
    for redshift_temp in redshift_space:
        sigma_space_mUV=[]
        for mUV_temp in mUV_space:
            fit=numpy.load('./COSMIC_VARIANCE_AT_HIGHZ/requirements/galaxy_cosmic_variance_fit_parameters_redshift_%.2f_stellar_mass_cut_%.2f_redshift_width_%.2f_aspect_ratio_1.00_log_no_of_points_7.npy'%(redshift_temp,mUV_temp,reference_redshift_width))
            sigma_space_mUV.append(10**linear(numpy.log10(length_in_Mpc),fit[0],fit[1]))           
        sigma_space_mUV=numpy.array(sigma_space_mUV)
        dat=scipy.interpolate.interp1d(mUV_space,numpy.log10(sigma_space_mUV),fill_value='extrapolate')
        sigma_space_redshift.append(10**dat(H_band_threshold))
    dat=scipy.interpolate.interp1d(redshift_space,numpy.log10(sigma_space_redshift),fill_value='extrapolate')
    
    fit=numpy.load('./COSMIC_VARIANCE_AT_HIGHZ/requirements/fit_params_redshift_width_dependence.npy')
    
    correction_for_redshift_width=10**linear(numpy.log10(redshift_width/reference_redshift_width),fit[0],fit[1])    
    return 10**dat(redshift)*correction_for_redshift_width


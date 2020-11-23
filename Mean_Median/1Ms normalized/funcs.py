import numpy as np
from scipy.stats import gamma
import scipy.integrate
import pandas
# TO DO: Update to most recent Planck results
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
from tqdm import tqdm

tot_sky = 41253.
file='extrapolate'
area_dict={'JADES Deep':46./3600/tot_sky, 
           'JADES Medium':190./3600/tot_sky, 
           'JADES Deep Extended':190./3600/tot_sky, 
           'JADES Medium Extended':5*190./3600/tot_sky,
           '1 pointing':  2*4.84/3600/tot_sky, 
           '20 pointings':  20*2*4.84/3600/tot_sky, 
           '100 pointings':  2*4.84/3600*100/tot_sky, 
           '1000 pointings':  2*4.84/3600*1000/tot_sky}

data = np.loadtxt("../erb2010_continuum.dat")
x = data[:,0]
y = data[:,1] * x**2 / 3.e18 # Convert to correct flux values

# Constants # Total degrees on the sky
masses = np.array([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0]) # Mass bins in log(M/Msun)

trials = 10000 # Number of trials for calculating median from gamma function
mmin = 9 # Minimum halo mass for HMF
mmax = 13.5 # Maximum halo mass for HMF

# Function to find closest wavelength in measured galaxy spectrum compared to given value
def nearest_wavelength(x, val):
    i = (np.abs(x - val)).argmin()
    return x[i]

# Define a class for different survey designs
class survey:
    def __init__(self, ABmax, area, cv_file, wavelength = 0, num_pointings = 1):
        self.ABmax = ABmax # Limiting AB magnitude of the survey
        self.area = area # Survey area
        self.file = cv_file # Name of cosmic variance table file
        self.point = num_pointings # Number of independent pointings
        self.wavelength = wavelength # Center wavelength for the chosen filter

# Method to take one trial of a gamma distribution with a given variance and mean
def trial(sig_v, mean, trials=10000):
    var = sig_v**2*mean**2
    k = mean**2/var
    t = var/mean
    rand = np.rint(np.random.gamma(k, scale=t, size=trials))
    return rand

#Song extrapolation

def lin(x,a=-0.04418885,b=-0.95909967):
    return a*x+b


if file=='extrapolate':
    def mean_median(survey, correction = True, path = "../../CV_Vals/"):
        # Get survey information
        import numpy as np
        from scipy.stats import gamma
        import scipy.integrate
        import pandas
        # TO DO: Update to most recent Planck results
        from astropy.cosmology import Planck18_arXiv_v2 as cosmo
        from tqdm import tqdm
        point = survey.point
        survey_area = survey.area * point
        ABmax = survey.ABmax
        # Find mean, median number of galaxies found for a given survey design
        # Absolute magnitude range over which to integrate
        # Range is larger than necessary to include all galaxies
        M = np.linspace(-30, 0, 1000)
        # Redshift values from CV calculations 
        df = pandas.read_csv(path+survey.file)
        z_vals = df["z"].values
        df.drop(columns=["z"])
        z_vals = z_vals[::-1]

        all_mean = []
        all_sig_v = []
        
        # Loop over z values
        for z in z_vals:
            # Correct magnitude to account for choice of telescope filter
            if(correction):
                # Get observed wavelength for a rest frame of 1600 Angstroms
                # Find closest wavelength to this value in the observed galaxy spectrum
                wl_obs = nearest_wavelength(x, 1600 * (1 + z))
                # Account for difference in observed magnitude due to difference in wavelength 
                # of the survey and the wavelength of the 
                ABmax_new = ABmax + 2.5 * np.log10(y[x==survey.wavelength]/y[x==wl_obs])
            else:
                ABmax_new = ABmax
            # Find volume of space within redshift bin
            dz = df.loc[df['z'] == z][["dz"]].values[0][0]
            vol = cosmo.comoving_volume(z + dz/2).value - cosmo.comoving_volume(z - dz/2).value

            # Apparent magnitude range over which to integrate
            m = M + cosmo.distmod(z).value

            # Schechter Extrapolation (Bouwens et al. 2015)
            Ms = -20.95 + 0.01 * (z - 6)
            p = 0.47 * 10**(-0.27 * (z - 6)) * 10**(-3)
            a = -1.87 - 0.1 * (z - 6)
            LF_center = np.log(10)/2.5 * p * (10**(0.4 * (Ms - M)))**(a+1) * np.exp(-10**(0.4 * (Ms - M)))

            means = []
            meds = []
            sigs = []
            # Loop over all mass bins
            for mass in masses:
                
                # Account for difference between Song and Bouwens wavelength
                diff = 2.5 * np.log10(y[x==nearest_wavelength(x,1505)]/y[x==nearest_wavelength(x,1605)])
                u_lim = 2. * (lin(z) - mass) + cosmo.distmod(z).value + diff # Dimmest object
                l_lim = 2. * (lin(z) - (mass + 0.5)) + cosmo.distmod(z).value + diff # Brightest object
                u_lim = min(u_lim, ABmax_new) # Compare dimmest object to telescope limit  
                # Apparent magnitude limits for given mass bin using Song et al. 2016
                # Find where apparent magnitude is within limits
                if mass == np.max(masses):
                    index = np.where(m <= u_lim)
                else:   
                    index = np.where((m >= l_lim) & (m <= u_lim))
                # Modify index for integration
                if index[0].size != 0:
                    index[0][-1] = index[0][-1] + 1

                # Integrate luminosity function over apparent magnitude to get number density
                int_center = np.trapz(LF_center[index], m[index])
                
                # Multiply by volume to get number
                # Get cosmic variance value for this z value and mass bin
                # CV values calculated as in Moster et al. 2010
                sig_v = float(df.loc[df['z'] == z][[str(mass)]].values[0][0]) / np.sqrt(int(point))
                sigs.append(sig_v)
                means.append(int_center * vol * survey_area / tot_sky)

            # Append mean and variance values of all masses for this redshift 
            all_mean.append(means)
            all_sig_v.append(sigs)  
        return all_mean, np.array(all_sig_v)

    #define sampling function
    def mm_sample(survey, correction = True, N=20, path = "../../CV_Vals/"):
        # Get survey information
        import numpy as np
        from scipy.stats import gamma
        import scipy.integrate
        import pandas
        from astropy.cosmology import Planck18_arXiv_v2 as cosmo
        from tqdm import tqdm
        point = survey.point
        survey_area = survey.area * point
        ABmax = survey.ABmax 
        df = pandas.read_csv(path+survey.file)
        # Absolute magnitude range over which to integrate
        # Range is larger than necessary to include all galaxies
        M = np.linspace(-30, 0, 1000)
        # Redshift values from CV calculations 
        z_vals = df["z"].values
        df.drop(columns=["z"])
        z_vals_orig = z_vals[::-1]
        avmean=np.zeros(len(z_vals))
        avmed=np.zeros(len(z_vals))
        for _ in tqdm(range(N)):
            z_vals = z_vals_orig*np.random.uniform(0.98,1.025,size=1)
            df["z"]=z_vals
            all_mean=[]
            all_sig_v=[]
            # Loop over z values
            for z in z_vals:
                # Correct magnitude to account for choice of telescope filter
                if(correction):
                    # Get observed wavelength for a rest frame of 1600 Angstroms
                    # Find closest wavelength to this value in the observed galaxy spectrum
                    wl_obs = nearest_wavelength(x, 1600 * (1 + z))
                    # Account for difference in observed magnitude due to difference in wavelength 
                    # of the survey and the wavelength of the 
                    ABmax_new = ABmax + 2.5 * np.log10(y[x==survey.wavelength]/y[x==wl_obs])
                else:
                    ABmax_new = ABmax
                
                # Find volume of space within redshift bin
                dz = df.loc[df['z'] == z][["dz"]].values[0][0]
                vol = cosmo.comoving_volume(z + dz/2).value - cosmo.comoving_volume(z - dz/2).value
                # print(vol)
                # Apparent magnitude range over which to integrate
                m = M + cosmo.distmod(z).value

                # Schechter Extrapolation (Bouwens et al. 2015)
                Ms = -20.95 + 0.01 * (z - 6)
                p = 0.47 * 10**(-0.27 * (z - 6)) * 10**(-3)
                a = -1.87 - 0.1 * (z - 6)
                LF_center = np.log(10)/2.5 * p * (10**(0.4 * (Ms - M)))**(a+1) * np.exp(-10**(0.4 * (Ms - M)))
                means = []
                meds = []
                sigs = []
                # Loop over all mass bins
                for mass in masses:
                    
                    # Account for difference between Song and Bouwens wavelength
                    diff = 2.5 * np.log10(y[x==nearest_wavelength(x,1505)]/y[x==nearest_wavelength(x,1605)])
                
                    # Apparent magnitude limits for given mass bin using Song et al. 2016
                    norm=lin(z)
                    u_lim = 2. * (norm - mass) + cosmo.distmod(z).value + diff # Dimmest object
                    l_lim = 2. * (norm - (mass + 0.5)) + cosmo.distmod(z).value + diff # Brightest object
                    u_lim = min(u_lim, ABmax_new) # Compare dimmest object to telescope limit  
                    
                    # Find where apparent magnitude is within limits
                    if mass == np.max(masses):
                        index = np.where(m <= u_lim)
                    else:   
                        index = np.where((m >= l_lim) & (m <= u_lim))

                    # Modify index for integration
                    if index[0].size != 0:
                        index[0][-1] = index[0][-1] + 1

                    # Integrate luminosity function over apparent magnitude to get number density
                    int_center = np.trapz(LF_center[index], m[index])
                    
                    # Multiply by volume to get number
                    means.append(int_center * vol * survey_area / tot_sky)
                    
                    # Get cosmic variance value for this z value and mass bin
                    # CV values calculated as in Moster et al. 2010
                    sig_v = float(df.loc[df['z'] == z][[str(mass)]].values[0][0]) / np.sqrt(int(point))
                    sigs.append(sig_v)

                # Append mean and variance values of all masses for this redshift 
                all_mean.append(means)
                all_sig_v.append(sigs)
                
            all_tot_mean = []
            all_med = []
            
            # Loop over all z values (high to low)
            for i in range(len(z_vals)):
                tot_mean = 0
                tot_trials = np.zeros(trials)
                # Find cumulative mean and median number of observed galaxies up to z
                for j in range(i+1):
                    means = all_mean[j]
                    sigs = all_sig_v[j]
                    # Run trials of gamma distribution to find median
                    for k in range(len(masses)):
                        if means[k] != 0:
                            one_trial = trial(sigs[k], means[k])
                            tot_trials = tot_trials + one_trial
                            tot_mean += means[k]
                all_tot_mean.append(tot_mean)
                med = np.median(tot_trials)
                all_med.append(med)  
            all_tot_mean, all_sig_v, all_med=np.array(all_mean), np.array(all_sig_v), np.array(all_med)
            avmean+=np.sum(all_tot_mean)
            avmed+=all_med
        return avmean/N, avmed/N, all_sig_v

    def mass_sample(survey, correction = True, steps=5, path = "../../CV_Vals/"):
            # Get survey information
            import numpy as np
            from scipy.stats import gamma
            import scipy.integrate
            import pandas
            # TO DO: Update to most recent Planck results
            from astropy.cosmology import Planck18_arXiv_v2 as cosmo
            from tqdm import tqdm
            point = survey.point
            survey_area = survey.area * point
            ABmax = survey.ABmax
            # Find mean, median number of galaxies found for a given survey design
            # Absolute magnitude range over which to integrate
            # Range is larger than necessary to include all galaxies
            M = np.linspace(-30, 0, 1000)
            # Redshift values from CV calculations 
            df = pandas.read_csv(path+survey.file)
            z_vals = df["z"].values
            z_vals = z_vals[::-1]

            all_mean = []
            all_sig_v = []
            all_indexlen = []
            # Loop over z values
            for z in z_vals:
                # Correct magnitude to account for choice of telescope filter
                if(correction):
                    # Get observed wavelength for a rest frame of 1600 Angstroms
                    # Find closest wavelength to this value in the observed galaxy spectrum
                    wl_obs = nearest_wavelength(x, 1600 * (1 + z))
                    # Account for difference in observed magnitude due to difference in wavelength 
                    # of the survey and the wavelength of the 
                    ABmax_new = ABmax + 2.5 * np.log10(y[x==survey.wavelength]/y[x==wl_obs])
                else:
                    ABmax_new = ABmax
                print(ABmax_new)
                # Find volume of space within redshift bin
                dz = df.loc[df['z'] == z][["dz"]].values[0][0]
                vol = cosmo.comoving_volume(z + dz/2).value - cosmo.comoving_volume(z - dz/2).value
                # print(vol)
                # Apparent magnitude range over which to integrate
                m = M + cosmo.distmod(z).value

                # Schechter Extrapolation (Bouwens et al. 2015)
                Ms = -20.95 + 0.01 * (z - 6)
                p = 0.47 * 10**(-0.27 * (z - 6)) * 10**(-3)
                a = -1.87 - 0.1 * (z - 6)
                LF_center = np.log(10)/2.5 * p * (10**(0.4 * (Ms - M)))**(a+1) * np.exp(-10**(0.4 * (Ms - M)))

                means = []
                sigs = []
                newmass = []
                indexlen = []
                # Loop over all mass bins
                for ma in masses:
                    sig1 = float(df.loc[df['z'] == z][[str(ma)]].values[0][0]) / np.sqrt(int(point))
                    if ma==11.0:
                        sig2=sig1
                    else:
                        sig2 = float(df.loc[df['z'] == z][[str(ma+0.5)]].values[0][0]) / np.sqrt(int(point))
                    for j in range(steps):
                        if ma==11.0:
                            mass=11.0
                        else:
                            mass=ma+0.5*(j/steps)
                        sig_v=sig1*(1-j/steps)+sig2*j/steps
                        # Account for difference between Song and Bouwens wavelength
                        # diff = 2.5 * np.log10(y[x==nearest_wavelength(x,1505)]/y[x==nearest_wavelength(x,1605)])
                        u_lim = 2. * (lin(z) - (mass-0.5)) + cosmo.distmod(z).value + diff # Dimmest object
                        l_lim = 2. * (lin(z) - (mass + 0.5)) + cosmo.distmod(z).value + diff # Brightest object
                        u_lim = min(u_lim, ABmax_new) # Compare dimmest object to telescope limit 
                        # print(z, mass, diff, u_lim, l_lim) 
                        # Apparent magnitude limits for given mass bin using Song et al. 2016
                        # Find where apparent magnitude is within limits
                        if mass == np.max(masses):
                            index = np.where(m <= u_lim)
                        else:   
                            index = np.where((m >= l_lim) & (m <= u_lim))
                        # Modify index for integration
                        if index[0].size != 0:
                            index[0][-1] = index[0][-1] + 1
                        indexlen.append(len(index[0]))
                        # Integrate luminosity function over apparent magnitude to get number density
                        int_center = np.trapz(LF_center[index], m[index])
                        print(z, mass, diff, u_lim, l_lim, len(index[0]), int_center)
                        # print(z, mass, int_center)
                        # print(mass, int_center)
                        # Multiply by volume to get number
                        # Get cosmic variance value for this z value and mass bin
                        # CV values calculated as in Moster et al. 2010
                        sigs.append(sig_v)
                        means.append(int_center * vol * survey_area / tot_sky)
                        newmass.append(mass)
                # Append mean and variance values of all masses for this redshift 
                all_mean.append(means[:-(steps)])
                all_sig_v.append(sigs[:-(steps)]) 
                all_indexlen.append(indexlen)
            all_mean=np.array(all_mean)
            all_mean[all_mean==0] = np.nan
            return np.array(all_mean), np.array(all_sig_v), newmass[:-steps], np.array(all_indexlen)

if file=='constant':
    def mean_median(survey, correction = True, path = "../../CV_Vals/"):
        # Get survey information
        point = survey.point
        survey_area = survey.area * point
        ABmax = survey.ABmax
        # Find mean, median number of galaxies found for a given survey design
        # Absolute magnitude range over which to integrate
        # Range is larger than necessary to include all galaxies
        M = np.linspace(-30, 0, 1000)
        # Redshift values from CV calculations 
        df = pandas.read_csv(path+survey.file)
        z_vals = df["z"].values
        df.drop(columns=["z"])
        z_vals = z_vals[::-1]

        all_mean = []
        all_sig_v = []
        
        # Loop over z values
        for z in z_vals:
            # Correct magnitude to account for choice of telescope filter
            if(correction):
                # Get observed wavelength for a rest frame of 1600 Angstroms
                # Find closest wavelength to this value in the observed galaxy spectrum
                wl_obs = nearest_wavelength(x, 1600 * (1 + z))
                # Account for difference in observed magnitude due to difference in wavelength 
                # of the survey and the wavelength of the 
                ABmax_new = ABmax + 2.5 * np.log10(y[x==survey.wavelength]/y[x==wl_obs])
            else:
                ABmax_new = ABmax
            # Find volume of space within redshift bin
            dz = df.loc[df['z'] == z][["dz"]].values[0][0]
            vol = cosmo.comoving_volume(z + dz/2).value - cosmo.comoving_volume(z - dz/2).value

            # Apparent magnitude range over which to integrate
            m = M + cosmo.distmod(z).value

            # Schechter Extrapolation (Bouwens et al. 2015)
            Ms = -20.95 + 0.01 * (z - 6)
            p = 0.47 * 10**(-0.27 * (z - 6)) * 10**(-3)
            a = -1.87 - 0.1 * (z - 6)
            LF_center = np.log(10)/2.5 * p * (10**(0.4 * (Ms - M)))**(a+1) * np.exp(-10**(0.4 * (Ms - M)))

            means = []
            sigs = []
            # Loop over all mass bins
            for mass in masses:
                # Account for difference between Song and Bouwens wavelength
                diff = 2.5 * np.log10(y[x==nearest_wavelength(x,1505)]/y[x==nearest_wavelength(x,1605)])
                if z<8:
                    u_lim = 2. * (-lin(z) - mass) + cosmo.distmod(z).value + diff # Dimmest object
                    l_lim = 2. * (-lin(z) - (mass + 0.5)) + cosmo.distmod(z).value + diff # Brightest object
                else:
                    u_lim = 2. * (-1.56 - mass) + cosmo.distmod(z).value + diff # Dimmest object
                    l_lim = 2. * (-1.56 - (mass + 0.5)) + cosmo.distmod(z).value + diff # Brightest object
                u_lim = min(u_lim, ABmax_new) # Compare dimmest object to telescope limit  
                # Apparent magnitude limits for given mass bin using Song et al. 2016
                # Find where apparent magnitude is within limits
                if mass == np.max(masses):
                    index = np.where(m <= u_lim)
                else:   
                    index = np.where((m >= l_lim) & (m <= u_lim))
                # Modify index for integration
                if index[0].size != 0:
                    index[0][-1] = index[0][-1] + 1

                # Integrate luminosity function over apparent magnitude to get number density
                int_center = np.trapz(LF_center[index], m[index])
                
                # Multiply by volume to get number
                # Get cosmic variance value for this z value and mass bin
                # CV values calculated as in Moster et al. 2010
                sig_v = float(df.loc[df['z'] == z][[str(mass)]].values[0][0]) / np.sqrt(int(point))
                sigs.append(sig_v)
                means.append(int_center * vol * survey_area / tot_sky)

            # Append mean and variance values of all masses for this redshift 
            all_mean.append(means)
            all_sig_v.append(sigs)  
        return all_mean, np.array(all_sig_v)

    def mm_sample(survey, correction = True, N=100, path = "../../CV_Vals/"):
        import numpy as np
        from scipy.stats import gamma
        import scipy.integrate
        import pandas
        from astropy.cosmology import Planck18_arXiv_v2 as cosmo
        from tqdm import tqdm
        # Get survey information
        point = survey.point
        survey_area = survey.area * point
        ABmax = survey.ABmax 
        df = pandas.read_csv(path+survey.file)
        # Absolute magnitude range over which to integrate
        # Range is larger than necessary to include all galaxies
        M = np.linspace(-30, 0, 1000)
        # Redshift values from CV calculations 
        z_vals = df["z"].values
        df.drop(columns=["z"])
        z_vals_orig = z_vals[::-1]
        avmean=np.zeros(len(z_vals))
        avmed=np.zeros(len(z_vals))
        for _ in tqdm(range(N)):
            z_vals = z_vals_orig*np.random.uniform(0.98,1.025,size=1)
            df["z"]=z_vals
            all_mean=[]
            all_sig_v=[]
            # Loop over z values
            for z in z_vals:
                # Correct magnitude to account for choice of telescope filter
                if(correction):
                    # Get observed wavelength for a rest frame of 1600 Angstroms
                    # Find closest wavelength to this value in the observed galaxy spectrum
                    wl_obs = nearest_wavelength(x, 1600 * (1 + z))
                    # Account for difference in observed magnitude due to difference in wavelength 
                    # of the survey and the wavelength of the 
                    ABmax_new = ABmax + 2.5 * np.log10(y[x==survey.wavelength]/y[x==wl_obs])
                else:
                    ABmax_new = ABmax
                
                # Find volume of space within redshift bin
                dz = df.loc[df['z'] == z][["dz"]].values[0][0]
                vol = cosmo.comoving_volume(z + dz/2).value - cosmo.comoving_volume(z - dz/2).value

                # Apparent magnitude range over which to integrate
                m = M + cosmo.distmod(z).value

                # Schechter Extrapolation (Bouwens et al. 2015)
                Ms = -20.95 + 0.01 * (z - 6)
                p = 0.47 * 10**(-0.27 * (z - 6)) * 10**(-3)
                a = -1.87 - 0.1 * (z - 6)
                LF_center = np.log(10)/2.5 * p * (10**(0.4 * (Ms - M)))**(a+1) * np.exp(-10**(0.4 * (Ms - M)))
                means = []
                meds = []
                sigs = []
                # Loop over all mass bins
                for mass in masses:
                    
                    # Account for difference between Song and Bouwens wavelength
                    diff = 2.5 * np.log10(y[x==nearest_wavelength(x,1505)]/y[x==nearest_wavelength(x,1605)])
                
                    # Apparent magnitude limits for given mass bin using Song et al. 2016
                    norm=lin(z)
                    u_lim = 2. * (-norm - mass) + cosmo.distmod(z).value + diff # Dimmest object
                    l_lim = 2. * (-norm - (mass + 0.5)) + cosmo.distmod(z).value + diff # Brightest object
                    u_lim = min(u_lim, ABmax_new) # Compare dimmest object to telescope limit  
                    
                    # Find where apparent magnitude is within limits
                    if mass == np.max(masses):
                        index = np.where(m <= u_lim)
                    else:   
                        index = np.where((m >= l_lim) & (m <= u_lim))

                    # Modify index for integration
                    if index[0].size != 0:
                        index[0][-1] = index[0][-1] + 1

                    # Integrate luminosity function over apparent magnitude to get number density
                    int_center = np.trapz(LF_center[index], m[index])
                    
                    # Multiply by volume to get number
                    means.append(int_center * vol * survey_area / tot_sky)
                    
                    # Get cosmic variance value for this z value and mass bin
                    # CV values calculated as in Moster et al. 2010
                    sig_v = float(df.loc[df['z'] == z][[str(mass)]].values[0][0]) / np.sqrt(int(point))
                    sigs.append(sig_v)

                # Append mean and variance values of all masses for this redshift 
                all_mean.append(means)
                all_sig_v.append(sigs)
                
            all_tot_mean = []
            all_med = []
            
            # Loop over all z values (high to low)
            for i in range(len(z_vals)):
                tot_mean = 0
                tot_trials = np.zeros(trials)
                # Find cumulative mean and median number of observed galaxies up to z
                for j in range(i+1):
                    means = all_mean[j]
                    sigs = all_sig_v[j]
                    # Run trials of gamma distribution to find median
                    for k in range(len(masses)):
                        if means[k] != 0:
                            one_trial = trial(sigs[k], means[k])
                            tot_trials = tot_trials + one_trial
                            tot_mean += means[k]
                all_tot_mean.append(tot_mean)
                med = np.median(tot_trials)
                all_med.append(med)  
            all_tot_mean, all_sig_v, all_med=np.array(all_mean), np.array(all_sig_v), np.array(all_med)
            avmean+=np.sum(all_tot_mean)
            avmed+=all_med
        return avmean/N, avmed/N, all_sig_v
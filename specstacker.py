from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import spectres
from astropy.coordinates import SkyCoord
import extinction as ext
from scipy import interpolate
import sfdmap as sfd
import spectres

from lmfit.models import GaussianModel
import warnings

# Ignore warnings
warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice")


line_list_dir='SDSS_emission_lines.txt'
dustmap_dir='dustmaps/sfddata-master'
filter_curve_dir='SDSS_filter_curves.txt'

def check_data(wave_spec, sigma_spec):
    ids=np.full(len(wave_spec), -99999)

    for i in range(len(wave_spec)):
        idx=np.where((wave_spec[i] >= 6600) & 
                          (wave_spec[i] <= 8380) & 
                          (sigma_spec[i] != np.inf))[0]
        if idx.size!=0:
            ids[i]=i
            
    if len(ids[ids==-99999])!=0:#i-band filter
        print('sample contains spectra with wavelengths below 6600A')
        print('number of ids to remove:', len(ids[ids==-99999]))
        return False, np.where(ids==-99999)[0]
    else:
        print('sample is checked')
        return True,ids


def correct_for_ext(wave_obs,flux,RA,DEC, dustmap_dir):
    """
   Correct for Galactic extinction.
   """
    # Get E(B-V) value from dustmap based on coordinates
    dustmap=sfd.SFDMap(dustmap_dir)#default scaling = 0.86 recalibration from Schlafly & Finkbeiner 2011
    E_BV=dustmap.ebv(RA,DEC)
        
    # Calculate Av value from E(B-V) and Rv
    Rv=3.1# the Rv value for the Fitzpatrick 1999 model
    Av=E_BV*Rv
    
    # Calculate extinction at each wavelength using Fitzpatrick 1999 model
    A_lambda=ext.fitzpatrick99(np.array(wave_obs,dtype='float'),a_v=Av,r_v=3.1,unit='aa')
    
    # Apply correction for Galactic extinction to the flux values
    flux_no_ext=flux*10**(A_lambda/2.5)
    
    # Return corrected flux values
    return flux_no_ext

def make_bins(wavs):
    """ 
    Given a series of wavelength points, find the edges and widths
    of corresponding wavelength bins. 
    """
    edges = np.zeros(wavs.shape[0]+1)
    widths = np.zeros(wavs.shape[0])
    edges[0] = wavs[0] - (wavs[1] - wavs[0])/2
    widths[-1] = (wavs[-1] - wavs[-2])
    edges[-1] = wavs[-1] + (wavs[-1] - wavs[-2])/2
    edges[1:-1] = (wavs[1:] + wavs[:-1])/2
    widths[:-1] = edges[1:-1] - edges[:-2]

    return edges, widths

def get_common_wave(wave_obs, z, sampling):
    """
    Define a common grid for resampling.
    """
    grid_rest = []  # Initialize empty list for grid arrays
    for i in range(len(wave_obs)):
        # Compute rest-frame wavelengths
        wave_rest = wave_obs[i] / (1 + z[i])
        # Compute grid for current rest-frame wavelength
        grid_i, _ = make_bins(wave_rest)
        grid_rest.append(grid_i)
    
    # Find minimum and maximum values in all grids
    grid_min, grid_max = np.inf, -np.inf
    for grid_i in grid_rest:
        grid_min = min(grid_min, np.min(grid_i))
        grid_max = max(grid_max, np.max(grid_i))
    
    # Define common grid as linearly spaced values within range
    grid_common=np.arange(np.ceil(grid_min),np.floor(grid_max+sampling),sampling)
    
    # Define wave_common as centers of common grid bins
    diff = grid_common[1] - grid_common[0]
    wave_common = np.arange(grid_common[0] + diff/2, grid_common[-1] + diff/2, diff)
    
    return wave_common


def resample(wave_common, wave_obs, flux, z):
    """
    Resample spectra onto a common grid.
    """
    # Initialize new flux array with NaNs
    flux_new = np.full((len(flux), len(wave_common)), np.nan)

    for i in range(len(flux)):
        # Only consider non-NaN values in flux
        no_nan = ~np.isnan(flux[i])

        # Compute rest-frame wavelengths
        wave_rest = wave_obs[i] / (1 + z[i])
        
        # Use spectres to resample flux onto common wavelength grid
        flux_new[i] = spectres.spectres(wave_common, wave_rest[no_nan], 
                                        flux[i][no_nan], fill=np.nan, verbose=False)

    return flux_new


def get_N(flux_new):
    """
    Compute the number distribution of the flux array.
    """
    # Create a copy of the flux array
    N = flux_new.copy()
    
    # Set all non-NaN values to 1.0
    N[~np.isnan(N)] = 1.0
    
    # Compute the sum of all values along the 0th axis
    N = np.nansum(N, axis=0)
    
    return N

def get_norm_range(wave, flux, z, sampling, line_list_dir=line_list_dir):
    """
    Get the normalized range of the spectrum, masking out emission lines.

    """
    # Convert the input wavelength array to a common wavelength array
    wave_common = get_common_wave(wave, z, sampling=1.0)

    # Resample the input flux array to the common wavelength array
    flux_new = resample(wave_common, wave, flux, z)

    # Compute the number distribution of the resampled array
    N = get_N(flux_new)
    
    # Define the range where all the spectra contribute to
    all_contr = np.where(N == max(N))[0]
    
    # Define the number of pixels to be used for the normalized range
    num_pix = int(np.round(len(wave_common) * 0.10))  # 10% of the data
    
    # Load the emission line wavelengths, names, and types from the input file
    em_list = np.loadtxt(line_list_dir, dtype='str', delimiter=',')
    line_wave = np.array([float(em_list[:,0][i]) for i in range(len(em_list))])
    line_name = np.array([str(em_list[:,1][i]) for i in range(len(em_list))])
    line_type = np.array([str(em_list[:,2][i]) for i in range(len(em_list))])
    
    # Mask out the emission lines from the normalized range
    norm_range = all_contr
    mask = np.full_like(norm_range, 0)
    
    for i in range(len(line_wave)):
        if 'broad' in line_type[i].lower():
            mask += ((wave_common[norm_range] < line_wave[i] - 100) |
                     (wave_common[norm_range] > line_wave[i] + 100))
        else:
            mask += ((wave_common[norm_range] < line_wave[i] - 50) |
                     (wave_common[norm_range] > line_wave[i] + 50))
        
    norm_range = norm_range[mask == len(line_wave)]
    
    # If the remaining normalized range is longer than the desired number of pixels,
    # truncate it to the last num_pix pixels
    if len(norm_range) >= num_pix:  
        norm_range = norm_range[len(norm_range) - num_pix:]
      
    # return the common wavelength array, the resamples and mormalised spectra, 
    # along with their number distribution and normalisation range    
    return wave_common, flux_new, N, norm_range

def Bootstrap(num,data):
    result=[]
    for i in range(num):
            sample=data[np.random.randint(low=0,high=len(data),size=len(data))]
            result.append(np.nanmedian(sample,axis=0))
    return np.nanstd(result,axis=0)


#ALGORITHM
def stacking_method(wave_obs, flux_obs, z, RA=None, DEC=None, 
                       line_list_dir =line_list_dir, sampling=1, dustmap_dir=dustmap_dir):
    """
    Stack a set of spectra onto a common grid, using a normalization range to 
    normalize each spectrum.
    """
    # Correct for Galactic extinction if RA and DEC are given
    if (RA is not None) and (DEC is not None):
        flux_no_ext = np.zeros_like(flux_obs)
        for i in range(len(flux_obs)):
            flux_no_ext[i] = correct_for_ext(wave_obs[i], flux_obs[i], RA[i], DEC[i], dustmap_dir)
        print('Spectra corrected for foreground extinction.')    
        flux_obs = flux_no_ext
      
    # Determine common wavelength grid, resample spectra onto common grid
    wave_common, flux_new, N, norm_range = get_norm_range(wave_obs, flux_obs, z, sampling, line_list_dir)
   
    # normalize each spectrum to a median value within the given range
    Fnorm = flux_new / np.median(flux_new[:, norm_range], axis=1, keepdims=True)
    
    # Stack spectra by taking the median of the normalized flux values
    F_stack = np.nanmedian(Fnorm, axis=0)
    
    # Bootstrap resampling to estimate uncertainties in the stacked spectrum
    F_Err = Bootstrap(1000, Fnorm)
    F_Err[N <= 5] = F_stack[N <= 5] / np.sqrt(N[N<=5])
    
    # Save stacked spectrum and associated information in a list
    stack = [F_stack, F_Err, N, wave_common, norm_range]
    
    return stack


def create_spectra(wave_spec, flux_spec, sigma_spec, z, flux_tem, wave_tem, filter_curve_dir):
    """
    Create spectra given a set of input parameters.
    """
    # Load i filter data
    filter_data = np.genfromtxt(filter_curve_dir, delimiter='', dtype='str', unpack=True, skip_header=4)
    idx_i = np.where(filter_data[4] != '...')[0]
    i_filter = np.array([float(filter_data[4][idx_i][i]) for i in range(len(filter_data[4][idx_i]))])
    i_wave = np.array([float(filter_data[0][idx_i][i]) for i in range(len(filter_data[0][idx_i]))])

    # Interpolate i filter data
    t_i = interpolate.interp1d(i_wave, i_filter)
    
    # Select only the spectra within the i-band filter's wavelength range and with finite sigma values
    idx_i = np.where((wave_spec >= i_wave[0]) & (wave_spec <= i_wave[-1]) & (sigma_spec != np.inf))[0]
    
    # Rescale template flux to match spectrum
    flux_new = spectres.spectres(wave_spec / (1 + z), wave_tem, flux_tem, spec_errs=None, fill=np.nan, verbose=False)
    
    # If there are spectra in the i-band filter range, calculate the i-band magnitude
    if idx_i.size != 0:
        F_i = sum(flux_spec[idx_i] * t_i(wave_spec[idx_i])) / sum(t_i(wave_spec[idx_i]))
        F_i_Jy = F_i * 10**-17 * (3.34 * 10**4.0 * 7480.0**2.0) # Convert to 10^-17 erg cm^-2 s^-1 A^-1
        mag_i = -2.5 * np.log10(F_i_Jy) + 8.9
    
        Fi_eff = sum(flux_new[idx_i] * t_i(wave_spec[idx_i])) / sum(t_i(wave_spec[idx_i]))
        F_i_Jy = 10 ** ((mag_i - 8.90) / -2.5)  # rearranging mAB_r=-2.5log10(F/Jy)+8.9
        F_i = F_i_Jy * 10 ** 17 / (3.34 * 10 ** 4.0 * 7480.0 ** 2.0)  # convert to SDSS units -> 10^-17 erg cm^-2 s^-1 A^-1 
        flux_denorm = flux_new * F_i / Fi_eff
    
        # Add noise to flux values
        flux_noise = flux_denorm + np.random.normal(0, 1, len(sigma_spec)) * sigma_spec
        flux_noise[sigma_spec == np.inf] = np.nan  # inf values can occur from ivar=0
    
        #return simulated flux
        return flux_noise
    
    else:
        print('spectrum not in i-band range')

def rescale(F,F_Err,wln,flux_tem,wave_tem,loc):
    """
    Rescale stack to match a template.
    """
    # resample the template to match the wavelength range of the input spectrum
    Ftem = spectres.spectres(wln, wave_tem, flux_tem, fill=np.nan, verbose=True)
    
    # rescale the flux and flux error of the input spectrum
    Fcal = F * np.average(Ftem[loc]) / np.average(F[loc])
    Fcal_Err = F_Err * np.average(Ftem[loc]) / np.average(F[loc])
    
    # return the rescaled flux and flux error
    return Fcal, Fcal_Err

def FD_method(x):
    """
    Computes the number of bins using the Freedman-Diaconis rule.
    """
    
    # Calculate the interquartile range
    q3, q1 = np.percentile(x, [75 ,25])
    iqr = q3 - q1
    
    # Calculate the bin width
    width=2*iqr*len(x)**(-1/3)
    
    # Return the number of bins
    return int((max(x)-min(x))/width)

def fit_gaussians_to_residual(stack_tem, stack_sim, loc):
    """
    Divides the residual obtained from two stacked spectra `stack_tem` and `stack_sim`. 
    For each division fits a Gaussian distribution and extracts the mean and standard deviation.
    """
    # Extract variables from input arrays
    fs, fs_std, N, wln = stack_sim[:4] # simulated spectrum variables
    wave_tem, ftem = stack_tem[3], stack_tem[0] # original/template spectrum variables
    
    # Rescale simulated spectrum to match the template
    fs, fs_std = rescale(fs, fs_std, wln, ftem, wave_tem, loc)
    
    # Compute chi-squared values between simulated and template spectra
    chi = (fs - ftem) / fs_std
    
    # Initialize variables for Gaussian fitting
    result = False
    parts = 30 # initial number of divisions of the wavelength range
    while not result:
        try:
            wave, mu, mu_err, sigma, sigma_err,  red_chi = [], [], [], [], [], [] # arrays to store fitted parameters
            
            #fix starting points
            wave.append(wln[0]) # initialize with first wavelength value
            mu.append(0) # initialize with mean of 0
            sigma.append(1) # initialize with stand. dev. of 1
            mu_err.append(0) # initialize with zero uncertainty
            sigma_err.append(0) # initialize with zero uncertainty
            red_chi.append(0) # initialize with value of 0
            
            num = (max(wln) - min(wln)) / parts # number of wavelength bins
            if num < 150: # if the number of bins is too small, raise an error
                raise ValueError
                
            # Loop over wavelength bins and fit Gaussian to each
            for i in range(parts):
                # Select wavelength range for current bin
                idx = np.where((wln >= wln[0] + i*num) & (wln < wln[0] + i*num + num))[0]
                
                # Compute histogram of chi-squared values in current range
                n, edge = np.histogram(chi[idx], bins=FD_method(chi[idx]))
                n = n.astype(float)
                centre = (edge[1:] + edge[:-1]) / 2
                centre = centre[n.nonzero()]
                n = n[n.nonzero()]
                
                # Fit Gaussian function to histogram
                gmodel = GaussianModel()
                params = gmodel.make_params(center=0, amplitude=max(n), sigma=1)
                result = gmodel.fit(n, params, x=centre, weights=1 / n**0.5)
                
                # If the fit has a high reduced chi-squared value, raise an error
                if result.redchi > 3.5:
                  #  print(i, 'rchi2=', result.redchi, len(idx))
                    raise ValueError
                
                # Store fitted parameters for current wavelength bin
                val = list(result.best_values.values())
                std = np.sqrt(np.diag(result.covar))
                mu.append(val[1])
                mu_err.append(std[1])
                sigma.append(val[2])
                sigma_err.append(std[2])
                wave.append(np.average(wln[idx]))
                red_chi.append(result.redchi)
                
            result = True # fitting was successful
                
        except ValueError:
            parts -= 1 # reduce number of wavelength bins and try again
            result = False
            
    print('Number of divisions used for the correction:', parts)      
    
    #update starting points
    mu_m = np.mean(mu[1:])
    sigma_m = np.mean(sigma[1:])
    mu[0] = mu_m
    sigma[0] = sigma_m
    
    #fix end points
    mu.append(mu_m)
    sigma.append(sigma_m)
    mu_err.append(0)
    sigma_err.append(0)
    wave.append(wln[-1])
    red_chi.append(0)
    
    return chi, wave, np.array(mu), np.array(mu_err), np.array(sigma), np.array(sigma_err), red_chi

def remove_nans(stack_tem, stack_sim):
    """
    Remove NaNs from the stacked spectra.
    """
    # Extract the data from the input stacks
    F1, F1_Err, N1, wln1, norm_range1 = stack_tem
    F2, F2_Err, N2, wln2, norm_range2 = stack_sim
    
    # Find the indices of the non-NaN values in both stacks
    idx = np.where((np.isnan(F1) == False) & (np.isnan(F2) == False))[0]
    
    # Adjust the normalization ranges based on the positions of the NaNs
    n_nans1 = np.sum(np.isnan(F1[:norm_range1[0]]))
    n_nans2 = np.sum(np.isnan(F2[:norm_range2[0]]))
    norm_range1 -= n_nans1
    norm_range2 -= n_nans2
    
    stack_tem = [F1[idx], F1_Err[idx], N1[idx], wln1[idx], norm_range1]
    stack_sim = [F2[idx], F2_Err[idx], N2[idx], wln2[idx], norm_range2]
    
    # Return the cleaned stacks
    return stack_tem, stack_sim


def correct_stack(stack,stack_sim):
    """
    Corrects the flux and associated error of a spectral stack based on a comparison with a simulated stack.
    """
    # remove any NaN values from the input stack
    stack, stack_sim = remove_nans(stack, stack_sim)
    
    # select the normalization range
    norm_range = stack[4]
    
    # fit the residual between the stack and the simulated stack with gaussians
    chi, wave, mu, mu_err, sigma, sigma_err, redchi = fit_gaussians_to_residual(stack, stack_sim, norm_range)
    
    # create interpolation functions for the fitted means and sigmas
    mu_func = interpolate.make_interp_spline(wave, mu)
    sigma_func = interpolate.make_interp_spline(wave, sigma)
    
    # correct the flux using the fitted mean and the simulated stack flux
    F_new = stack[0] + mu_func(stack[3]) * stack_sim[1]
    
    # correct the flux uncertainty using the fitted sigma and the simulated stack flux uncertainty
    F_Err_new = stack[1] * sigma_func(stack[3])
    
    # create a new stack with the corrected flux and flux uncertainty
    stack_new = list(stack).copy()
    stack_new[0] = F_new
    stack_new[1] = F_Err_new
    
    # return the corrected stack
    return stack_new


def stack_spectra(wave_spec,flux_spec,sigma_spec,z,RA,DEC,zbins=None,sampling=1, 
                    line_list_dir=line_list_dir, dustmap_dir=dustmap_dir,
                    filter_curve_dir=filter_curve_dir):#correction='yes'
    """
    Calculates corrected stacks for given redshift bins.

    Parameters:
        -----------
    wave_spec : array_like 
        1D array of wavelengths for all spectra.
    flux_spec : array_like 
        2D array of fluxes for all spectra.
    sigma_spec : array_like 
        2D array of flux errors for all spectra.
    z : array_like
        1D array of redshifts for all spectra.
    RA : array_like
        1D array of right ascension coordinates for all spectra, if foreground correction needed.
    DEC : array_like
       1D array of declination coordinates for all spectra, if foreground correction needed.
    zbins : list
        List of tuples where each tuple contains the lower and upper bound of the redshift bin.
    sampling : int
        Number of samples to use for resampling.
    line_list_dir :
        Directory of a list of spectral lines that includes 
        the line central wavelegths, the line name and the line type  
    dustmap_dir : 
        Directory of the dustmaps  
    filter_curve_dir : 
        Directory of the filter curve used to simulate the spectra  
    

    Returns:
    --------
    stacks : list
        Dictionary of original stacks.
    stacks_sim : list
        Dictionary of simulated stacks.
    staacks_corr : list
        Dictionary of corrected stacks.
    spec_idx : list
        Dictionary of indexes for spectra in each redshift bin.
    """
    #check_data
    # check, ids = check_data(wave_spec, sigma_spec)
    # if check==True or check==False:
    
    #initialize empty lists to store the indexes, stacks, simulated stacks, and corrected stacks
    #idx, stacks, stacks_sim, corrs = [], [], [], [] 
    idx=[]
    spec_idx, stacks, stacks_sim, stacks_corr={},{},{},{}
        
    
    # if zbins not specified take the min amx max z value
    if zbins is None:
            zmin, zmax = np.min(z), np.max(z)
            zbins = [(zmin, zmax)]
    print('zbin:', zbins)
    for i in range(len(zbins)):
            #select indexes of spectra within the redshift bin and print the total selected
            idx.append(np.where((z >= zbins[i][0]) & (z <= zbins[i][1]))[0]) 
            print('Number of spectra to stack = ', len(idx[i]))
        
            spec_idx['zbin='+str(i)] = {'spec_idx': idx[i]}
            # select wavelegth, flux, flux error, redshift, RA, DEC based on the selected indexes
            wave, flux, sigma = wave_spec[idx[i]], flux_spec[idx[i]], sigma_spec[idx[i]] 
            z_i = z[idx[i]] 
            RA_i = RA[idx[i]]
            DEC_i = DEC[idx[i]] 
        
            # stack the spectra
            stack = stacking_method(wave, flux, z_i, RA=RA_i, DEC=DEC_i, 
                                       sampling=sampling, line_list_dir=line_list_dir, dustmap_dir=dustmap_dir) 
            #  stacks.append(stack) 
            stacks['zbin='+str(i)] = {'flux': stack[0], 'flux_err': stack[1], 'N': stack[2],
                                  'wln': stack[3], 'norm_range': stack[4]}

            # create simulated spectra based on the stacking result
            spec=[]
            for j in range(len(idx[i])):
                    spec.append(create_spectra(wave[j], flux[j], sigma[j],
                                              z_i[j], stack[0], stack[3], filter_curve_dir))
            # stack the simulated spectra
            stack_sim = stacking_method(wave, spec, z_i, RA=None, DEC=None, 
                                       sampling=sampling, line_list_dir=line_list_dir, dustmap_dir=dustmap_dir)
            stacks_sim['zbin='+str(i)] = {'flux': stack_sim[0], 'flux_err': stack_sim[1], 'N': stack_sim[2],
                                      'wln': stack_sim[3], 'norm_range': stack_sim[4]}
    
            #correct the original stack using the simulated stack
            corr=correct_stack(stack, stack_sim)
            stacks_corr['zbin='+str(i)] = {'flux': corr[0], 'flux_err': corr[1], 'N': corr[2],
                                      'wln': corr[3], 'norm_range': corr[4]}
        
        #return the original stacks, simulated stacks, corrected stacks, and the indexes of selected spectra
    return stacks, stacks_sim, stacks_corr, spec_idx 

 

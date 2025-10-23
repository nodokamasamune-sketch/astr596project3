
def read_draine_opacity(filename):
    """
    Read Draine dust opacity file.
    
    Returns:
    --------
    data : dict
        'wavelength': array in cm
        'kappa': array in cm^2/g
        'albedo': array
        'dust_to_H': float (mass ratio)
    """
    import numpy as np
    dust_data = np.genfromtxt(filename, skip_header=48, names=['wavelength', 'albedo', 'g', 'C_ext_H', 'K_abs'])
    
    dust_per_H = 1.798e-26
    
    data = {}
    data['wavelength'] = dust_data['wavelength'] * (10**(-4)) #CONVERT MICRONS TO CM
    data['kappa'] = dust_data['C_ext_H'] / dust_per_H #cm^2 / gram
    data['albedo'] = dust_data['albedo']
    data['dust_to_H'] = 0.01 #percent ?? --> should it be 1?
    
    return data

def calculate_band_averaged_opacity(draine_data, band_min, band_max, t):
    """
    Calculate Planck mean opacity for given stellar temperature.
    
    Parameters:
    -----------
    draine_data : dict
        From read_draine_opacity
    band_min, band_max : float
        Band limits in cm
    T_star : float
        Stellar temperature (K)
    
    Returns:
    --------
    kappa_avg : float
        Band-averaged opacity (cm^2/g)
    """
    import numpy as np
    from scipy import integrate
    from utils import planck_function
    import constants as cons
    
    lambda_ = np.linspace(band_min, band_max, 1000)
    

    kappa_interp = np.interp(lambda_, draine_data['wavelength'], draine_data['kappa'])
    b_lambda = planck_function(lambda_, t)
    kappa_int = integrate.simpson(kappa_interp * b_lambda, lambda_)
    b_int = integrate.simpson(b_lambda, lambda_)
    kappa_avg = kappa_int / b_int
    
    return kappa_avg


class Band:
    def __init__(self, name, lambda_min, lambda_max):
        """
        Wavelength band definition.
        
        Required attributes:
        - name: 'B', 'V', or 'K'
        - lambda_min: minimum wavelength (cm)
        - lambda_max: maximum wavelength (cm)
        - lambda_center: central wavelength (cm)
        """
        pass
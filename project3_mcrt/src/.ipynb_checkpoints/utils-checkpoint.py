def planck_function(wavelength, temperature):
    """
    Planck function B_lambda(T) in CGS units.
    
    Parameters:
    -----------
    wavelength : float or array
        Wavelength in cm
    temperature : float
        Temperature in K
    
    Returns:
    --------
    B_lambda : float or array
        Planck function in erg/s/cm^2/sr/cm
    """
    import numpy as np
    import constants as cons

    B_lambda = ((2 * cons.h * cons.CSOL**2) / wavelength**5) * 1/(np.exp((cons.h * cons.CSOL)/ (wavelength * cons.k * temperature)) -1)
    return B_lambda

def integrate_band(func, band_min, band_max, *args):
    """
    Integrate function over wavelength band. Use scipy.integrate.simpson.
    
    Parameters:
    -----------
    func : callable
        Function to integrate
    band_min, band_max : float
        Band limits in cm
    *args : additional arguments for func
    
    Returns:
    --------
    result : float
        Integrated value
    """

    from scipy import integrate
    import numpy as np

    lambda_ = np.linspace(band_min, band_max, 1000)
    B = func(lambda_, *args)
    int_B = integrate.simpson(B, lambda_)                         

    return int_B
# zams.py
"""
ZAMS mass-luminosity and mass-radius relations from Tout et al. (1996)
Self-contained functions with coefficients from the paper
Valid for masses 0.1 - 100 M_sun and Z = 0.0001 - 0.03
"""
import numpy as np

def luminosity(M, Z=0.02, solar_Z_only=True):
    """
    Calculate ZAMS luminosity using Tout et al. (1996) Eq. 1
    
    L/L_sun = (α M^2.5 + β M^11) / (M^3 + γ + δ M^5 + ε M^7 + ζ M^8 + η M^9.5)
    
    Parameters
    ----------
    M : float or np.ndarray
        Stellar mass in solar masses
        Valid range: 0.1 - 100 M_sun
    Z : float, optional
        Metallicity (default: 0.02 for solar)
        Valid range from Tout et al.: 0.0001 - 0.03
    solar_Z_only : bool, optional
        If True (default), only Z=0.02 is implemented
        Set to False for extension with Z-dependence
    
    Returns
    -------
    L : float or np.ndarray
        Luminosity in solar luminosities
    
    Raises
    ------
    AssertionError
        If mass or metallicity outside valid range
    
    References
    ----------
    Tout et al. (1996) MNRAS 281, 257
    See equations (1) for luminosity formula
    See equations (3)-(4) for metallicity dependence
    """
    # Input validation for mass
    # isinstance(M, np.ndarray) checks if M is a numpy array object
    # This lets our function handle both single values and arrays!
    if isinstance(M, np.ndarray):
        # np.all() returns True only if ALL elements satisfy the condition
        assert np.all(M >= 0.1) and np.all(M <= 100), \
            "Mass must be between 0.1 and 100 M_sun (Tout et al. 1996 validity range)"
    else:
        assert 0.1 <= M <= 100, \
            "Mass must be between 0.1 and 100 M_sun (Tout et al. 1996 validity range)"
    
    # Assert Z is in Tout et al. valid range (0.0001 to 0.03)
    if isinstance(Z, np.ndarray):
        # np.all() returns True only if ALL elements satisfy the condition
        assert np.all(Z >= 0.0001) and np.all(Z <= 0.03), \
            "Metallicity must be between 0.0001 and 0.03 (Tout et al. 1996 validity range)"
    else:
        assert 0.0001 <= Z <= 0.03, \
            "Metallicity must be between 0.0001 and 0.03 (Tout et al. 1996 validity range)"
    # If solar_Z_only is True, also assert Z == 0.02
    if solar_Z_only==True:
        assert Z == 0.02, \
            'Z must equal 0.02.'
        
    # Coefficients from Table 1 (check equations 3-4 to understand the table structure)
    a = 0.39704170
    b = 8.527626
    g = 0.00025546
    d = 5.4328890
    e = 5.563569
    z = 0.7886606
    et = 0.00586685
    
    # Equation (1) from Tout et al. (1996) to calculate luminosity
    L = (a*(M**5.5) + b*(M**11)) / (g + M**3 + d*(M**5) + e*(M**7) + z*(M**8) + et*(M**9.5))
    
    # IMPORTANT: Note the fractional exponents (M^2.5, M^9.5, etc.)!
    return L

# Implement radius() function following the same pattern
# Use Equation (2) and coefficients from Table 2
# Check equations (3)-(4) to understand which column corresponds to Z=0.02
def radius(M, Z=0.02, solar_Z_only=True):
    if isinstance(M, np.ndarray):
        # np.all() returns True only if ALL elements satisfy the condition
        assert np.all(M >= 0.1) and np.all(M <= 100), \
            "Mass must be between 0.1 and 100 M_sun (Tout et al. 1996 validity range)"
    else:
        assert 0.1 <= M <= 100, \
            "Mass must be between 0.1 and 100 M_sun (Tout et al. 1996 validity range)"
    
    # Assert Z is in Tout et al. valid range (0.0001 to 0.03)
    if isinstance(Z, np.ndarray):
        # np.all() returns True only if ALL elements satisfy the condition
        assert np.all(Z >= 0.0001) and np.all(Z <= 0.03), \
            "Metallicity must be between 0.0001 and 0.03 (Tout et al. 1996 validity range)"
    else:
        assert 0.0001 <= Z <= 0.03, \
            "Metallicity must be between 0.0001 and 0.03 (Tout et al. 1996 validity range)"
    # If solar_Z_only is True, also assert Z == 0.02
    if solar_Z_only==True:
        assert Z == 0.02, \
            'Metallicity must equal 0.02.'

        
    # Coefficients from Table 2 (check equations 3-4 to understand the table structure)
    t = 1.715359
    i = 6.597788
    k = 10.08855
    l = 1.012495
    m = 0.07490166
    v = 0.01077422
    e = 3.082234
    o = 17.84778
    p = 0.00022582

    # Equation (2) from Tout et al. (1996) to calculate radius
    R = (t*(M**2.5) + i*(M**6.5) + k*(M**11) + l*(M**19) + m*(M**19.5)) / (v + e*(M**2) + o*(M**8.5) + M**18.5 + p*(M**19.5))

    return R


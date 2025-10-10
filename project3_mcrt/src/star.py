
class Star:
    def __init__(self, mass, Z=0.02):
        """
        Initialize ZAMS star.
        
        Required attributes:
        - mass: stellar mass in M_sun
        - Z: metallicity
        - T_eff: effective temperature (K)
        - L_bol: bolometric luminosity (erg/s)
        - R: stellar radius (cm)
        - x, y, z: position (cm) NEW
        - id: unique identifier NEW
        - kappa_band: dict of band-averaged opacities {'B': value, 'V': value, 'K': value} NEW
        - L_band: dict of band luminosities {'B': value, 'V': value, 'K': value} NEW
        
        IMPORTANT: Calculate opacities and luminosities HERE during initialization!
        After computing T_eff from ZAMS:
            self.kappa_band = {}
            self.L_band = {}
            for band in ['B', 'V', 'K']:
                self.kappa_band[band] = self.calculate_planck_mean_opacity(...)
                self.L_band[band] = self.calculate_band_luminosity(...)
        """
        # Use your Project 1 ZAMS calculations
        import numpy as np
        from zams import luminosity, radius
        import constants as cons

        self.mass = mass
        assert 0.1 <= self.mass <= 100, \
        "the mass must be between 0.1 and 100 solar masses."
        
        self.luminosity = luminosity(mass)
        self.radius = radius(mass)

        self.t_eff = ((self.luminosity)**(0.25) * self.radius**(-0.5)) * cons.T_sun        

        
    
    def calculate_band_luminosity(self, band_min, band_max):
        """
        Calculate luminosity in specific band.
        
        Returns:
        --------
        L_band : float
            Luminosity in band (erg/s)
        """
        from utils import planck_function, integrate_band
        import constants as cons
        import numpy as np

        luminosity = self.luminosity
        t = self.t_eff
        int_L = integrate_band(planck_function, band_min, band_max, t)
        L_b = luminosity * int_L / ((cons.SIGMA_SB * t**4) / np.pi)
        
        return L_b
    
    def calculate_planck_mean_opacity(self, draine_data, band_min, band_max):
        """
        Calculate Planck mean opacity for this star's temperature.
        
        Returns:
        --------
        kappa : float
            Planck mean opacity (cm^2/g)
        """
        pass
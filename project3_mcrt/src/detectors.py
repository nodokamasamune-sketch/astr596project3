# Students: Add your own import statements here
# Suggested: numpy

class EscapeTracker:
    def __init__(self, bands):
        """
        Track escaping packets.
        
        Required attributes:
        - L_escaped_total: total escaped luminosity
        - L_escaped_by_band: dict of escaped luminosity per band
        - n_escaped_by_band: dict of packet counts per band
        - escape_directions: list of (theta, phi) for each packet
        """
        pass
    
    def record_escape(self, packet):
        """
        Record an escaping packet.
        """
        pass
    
    def calculate_escape_fractions(self, L_input_by_band):
        """
        Calculate escape fraction for each band.
        
        Returns:
        --------
        f_escape : dict
            Escape fraction for each band
        """
        pass

def create_escape_map(escape_directions, n_theta=180, n_phi=360):
    """
    Create 2D histogram of escape directions.
    
    Returns:
    --------
    escape_map : 2D array
        Map of escaping luminosity vs direction
    """
    pass
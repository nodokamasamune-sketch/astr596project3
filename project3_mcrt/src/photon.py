
class Photon:
    def __init__(self):
        """
        Initialize photon packet.
        
        Required attributes:
        - x, y, z: position (cm)
        - dir_x, dir_y, dir_z: direction unit vector
        - L: luminosity carried (erg/s)
        - band: band identifier ('B', 'V', or 'K')
        - star_id: source star ID
        - kappa: opacity for this packet (cm^2/g)
        """
        pass
    
    def move(self, distance):
        """
        Move photon along its direction.
        
        Parameters:
        -----------
        distance : float
            Distance to move (cm)
        """
        pass

def emit_packet_from_star(star, band, L_packet):
    """
    Create packet emitted from stellar surface.
    
    Parameters:
    -----------
    star : Star object
    band : Band object or string
    L_packet : float
        Luminosity carried by packet
    
    Returns:
    --------
    packet : Photon object
    """
    pass

def sample_isotropic_direction():
    """
    Sample random isotropic direction.
    
    Returns:
    --------
    dir_x, dir_y, dir_z : float
        Unit direction vector
    """
    pass

def distance_to_next_boundary(packet, grid):
    """
    Calculate distance to next cell boundary.
    
    Parameters:
    -----------
    packet : Photon object
    grid : Grid object
    
    Returns:
    --------
    d_next : float
        Distance to next boundary (cm)
    face : str
        Which boundary ('x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max')
    """
    pass

def propagate_packet(packet, grid):
    """
    Propagate packet through grid until absorbed or escaped.
    
    Parameters:
    -----------
    packet : Photon object
    grid : Grid object
    
    Returns:
    --------
    outcome : str
        'absorbed' or 'escaped'
    location : tuple
        (ix, iy, iz) if absorbed, (x, y, z) if escaped
    """
    pass
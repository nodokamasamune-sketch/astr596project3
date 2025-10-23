
class Photon:
    def __init__(self, star, band):
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
        
        # randomly generated position on stellar surface

        a = np.random.uniform(0, 1)
        b = np.random.uniform(0, 1)

        theta = np.arccos(2*a - 1)
        phi = 2 * np.pi * b

        self.x = star.x + star.radius * np.sin(theta) * np.cos(phi)
        self.y = star.y + star.radius * np.sin(theta) * np.sin(phi)
        self.z = star.z + star.radius * np.cos(theta)

        # randomly generated direction
        self.x_dir, self.y_dir, self.z_dir = sample_isotropic_direction()
        
        
        # luminosity carried
        N = 1 # number of photon packets
        self.L = 

        # band
        self.band = band

        # source star ID
        self.star_id = star.id

        #kappa
        self.kappa = star.kappa_band[band]

        
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

    photon = Photon(star, band)

    return photon

    

def sample_isotropic_direction():
    """
    Sample random isotropic direction.
    
    Returns:
    --------
    dir_x, dir_y, dir_z : float
        Unit direction vector
    """
    # randomly generated direction unit vector

    c = np.random.uniform(0, 1)
    d = np.random.uniform(0, 1)
    theta_dir = np.arccos(2*c - 1)
    phi_dir = 2 * np.pi * d
        
    x_dir = np.sin(theta_dir) * np.cos(phi_dir)
    y_dir = np.sin(theta_dir) * np.sin(phi_dir)
    z_dir = np.cos(theta_dir)

    return x_dir, y_dir, z_dir
    

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

    x = packet.x
    y = packet.y
    z = packet.z

    x_dir = packet.x_dir
    y_dir = packet.y_dir
    z_dir = packet.z_dir
    
    ix, iy, iz = get_cell_indices(self, x, y, z)

    # Compute distance to next x-boundary
    if x_dir > 0:
        x_edge = grid.x_min + (ix + 1) * grid.dx
        dx = (x_edge - x) / x_dir
        face_x = 'x_max'
    elif x_dir < 0:
        x_edge = grid.x_min + ix * grid.dx
        dx = (x_edge - x) / x_dir
        face_x = 'x_min'

    # Compute distance to next y-boundary
    if y_dir > 0:
        y_edge = grid.y_min + (iy + 1) * grid.dy
        dy = (y_edge - y) / y_dir
        face_y = 'y_max'
    elif y_dir < 0:
        y_edge = grid.y_min + iy * grid.dy
        dy = (y_edge - y) / y_dir
        face_y = 'y_min'


    # Compute distance to next z-boundary
    if z_dir > 0:
        z_edge = grid.z_min + (iz + 1) * grid.dz
        dz = (z_edge - z) / z_dir
        face_z = 'z_max'
    elif z_dir < 0:
        z_edge = grid.z_min + iz * grid.dz
        dz = (z_edge - z) / z_dir
        face_z = 'z_min'


    # Determine minimum distance
    distances = [dx, dy, dz]
    faces = [face_x, face_y, face_z]
    min_index = distances.index(min(distances))

    return distances[min_index], faces[min_index]


    


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

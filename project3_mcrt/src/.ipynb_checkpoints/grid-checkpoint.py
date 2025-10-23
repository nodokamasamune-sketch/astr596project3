# Students: Add your own import statements here  
# Suggested: numpy
import numpy as np


class Grid:
    def __init__(self, n_cells=128, L_pc=1.0, f_dust_to_gas=0.01):
        """
        Initialize 3D grid.
        
        Parameters:
        -----------
        n_cells : int
            Number of cells per dimension (nx = ny = nz)
        L_pc : float
            Box size in parsecs
        f_dust_to_gas : float
            Dust-to-gas mass ratio (default 0.01)
        
        Required attributes:
        - nx, ny, nz: number of cells
        - L: box size in cm
        - dx, dy, dz: cell size in cm
        - x_min, x_max: boundaries (-L/2, L/2)
        - y_min, y_max: boundaries (-L/2, L/2)
        - z_min, z_max: boundaries (-L/2, L/2)
        - f_dust_to_gas: dust-to-gas ratio
        - rho_gas: 3D array of gas density
        - L_absorbed: 3D array for absorbed luminosity
        - n_absorbed: 3D array for packet count
        """
        self.nx = n_cells
        self.ny = n_cells
        self.nz = n_cells

        
        self.L = L_pc * 3.086e18 #convert pc to cm (to maintain CGS units)

        self.dx = self.L / self.nz
        self.dy = self.L / self.ny
        self.dz = self.L / self.nz
        
        self.min = -self.L / 2
        self.max = self.L / 2

        self.x_min = self.min
        self.x_max = self.max
        self.y_min = self.min
        self.y_max = self.max
        self.z_min = self.min
        self.z_max = self.max

        self.xvalues = np.linspace(self.x_min, self.x_max, self.nx)
        self.yvalues = np.linspace(self.y_min, self.y_max, self.ny)
        self.zvalues = np.linspace(self.z_min, self.z_max, self.nz)
        
        self.X, self.Y, self.Z = np.meshgrid(self.xvalues, self.yvalues, self.zvalues, indexing='ij')
        self.grid = np.stack((self.X, self.Y, self.Z), axis=-1)

        
    
    def get_cell_indices(self, x, y, z):
        """
        Get cell indices for position.
        
        Returns:
        --------
        ix, iy, iz : int
            Cell indices
        """
        ix = int((x - self.x_min) / self.dx)
        iy = int((y - self.y_min) / self.dy)
        iz = int((z - self.z_min) / self.dz)

        #consider adding check to make sure ix iy iz are within bounds
        
        return ix, iy, iz
        
    
    def is_inside(self, x, y, z):
        """
        Check if position is inside grid.
        
        Returns:
        --------
        inside : bool
        """
        
        inside = True
        if x < self.min or x > self.max:
            inside = False
        if y < self.min or y > self.max:
            inside = False
        if z < self.min or z > self.max:
            inside = False
        
        return inside
                
    
    def get_dust_density(self, ix, iy, iz):
        """
        Get dust density in cell.
        
        Returns:
        --------
        rho_dust : float
            Dust density (g/cm^3)
        """
        #dust density uniform throughout grid
        
        rho_dust = 3.84e-23 #g / cm^3

        return rho_dust
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

        self.min = -self.L / 2
        self.max = self.L / 2

        self.xvalues = np.linspace(self.min, self.max, self.nx)
        self.yvalues = np.linspace(self.min, self.max, self.ny)
        self.zvalues = np.linspace(self.min, self.max, self.nz)
        
        self.X, self.Y, self.Z = np.meshgrid(self.xvalues, self.yvalues, self.zvalues)
        self.grid = np.stack(self.X, self.Y, self.Z)
        pass
    
    def get_cell_indices(self, x, y, z):
        """
        Get cell indices for position.
        
        Returns:
        --------
        ix, iy, iz : int
            Cell indices
        """
        
        ix = np.where(self.X == x)
        iy = np.where(self.Y == y)
        iz = np.where(self.Z == z)
        
        
        
        
        pass
    
    def is_inside(self, x, y, z):
        """
        Check if position is inside grid.
        
        Returns:
        --------
        inside : bool
        """
        pass
    
    def get_dust_density(self, ix, iy, iz):
        """
        Get dust density in cell.
        
        Returns:
        --------
        rho_dust : float
            Dust density (g/cm^3)
        """
        pass
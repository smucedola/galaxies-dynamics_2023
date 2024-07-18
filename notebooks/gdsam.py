import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from tqdm.notebook import tqdm
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.stats import gaussian_kde




def distance(x1, y1, z1,
                 x2, y2, z2):
        
        return np.sqrt((x1-x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)


def center_of_mass(x, y, z, m):
    M = np.sum(m)
    x_cm = np.sum(x * m[:, None], axis=0) / M
    y_cm = np.sum(y * m[:, None], axis=0) / M
    z_cm = np.sum(z * m[:, None], axis=0) / M
    
    return x_cm, y_cm, z_cm


class initial_data:
    """
    Initialize data to standard distributions.
    """

    def __init__(self):
        """
        Constructor of sh1_init.
        """
        # Initialize attributes to None
        self.N = None
        self.m = None
        self.r, self.theta, self.phi, self.z = None, None, None, None
        self.x, self.y, self.z = None, None, None
        self.vx, self.vy, self.vz = None, None, None

    def hom_sph(self, N, M_tot, radius):
        """
        Create initial positions and velocities of a homogeneous spherical distribution.

        Args:
            N (int): Number of particles.
            M_tot (float): Total mass.
            radius (float): Radius of the sphere.
            v (tuple): Initial velocities. Default is (0, 0, 0).
        """
        self.N = N
        self.m = np.full(N, M_tot / N)
        
        # Define lambda functions for inverse CDFs
        I_r = lambda u: u**(1/3) * radius
        I_theta = lambda u: np.arccos(2 * u - 1)
        
        # Sample random numbers for r and theta
        u_r = np.random.uniform(0, 1, N)
        u_theta = np.random.uniform(0, 1, N)
        
        # Generate positions following homogeneous spherical distribution
        self.r = I_r(u_r)
        self.theta = I_theta(u_theta)
        self.phi = np.random.uniform(0, 2 * np.pi, N)
        
        self.x = self.r * np.sin(self.theta) * np.cos(self.phi)
        self.y = self.r * np.sin(self.theta) * np.sin(self.phi)
        self.z = self.r * np.cos(self.theta)


    def plummer_model(self, N, M_tot, radius):
        """
        Create initial positions and velocities of a Plummer model distribution.

        Args:
            N (int): Number of particles.
            M_tot (float): Total mass.
            radius (float): Plummer radius.
        """
        self.N = N
        self.m = np.full(N, M_tot / N)
        # Define lambda functions for inverse CDFs
        I_r = lambda u: radius* np.sqrt(u**(2/3)/(1- u**(2/3)))
        I_theta = lambda u: np.arccos(2 * u - 1)
        
        # Sample random numbers for r and theta
        u_r = np.random.uniform(0, 1, N)
        u_theta = np.random.uniform(0, 1, N)
        
        # Generate positions following Plummer model distribution
        self.r = I_r(u_r)
        self.theta = I_theta(u_theta)
        self.phi = np.random.uniform(0, 2 * np.pi, N)
        
        self.x = self.r * np.sin(self.theta) * np.cos(self.phi)
        self.y = self.r * np.sin(self.theta) * np.sin(self.phi)
        self.z = self.r * np.cos(self.theta)

    
    def disk(self, N, M_tot, radius):
        self.N = N
        self.m = np.full(N, M_tot / N)
        # Define lambda functions for inverse CDFs
        I_r = lambda u: radius*np.sqrt((1/(1-u)**2-1))
        I_theta = lambda u: 2*np.pi*u
        
        # Sample random numbers for r and theta
        u_r = np.random.uniform(0, 1, N)
        u_theta = np.random.uniform(0, 1, N)
        
        # Generate positions following Plummer model distribution
        self.r = I_r(u_r)
        self.theta = I_theta(u_theta)-np.pi
        self.z = np.zeros(self.N) #this is now z in rality, cyl coords

        self.x = self.r * np.cos(self.theta)
        self.y = self.r * np.sin(self.theta)
        self.z = self.z


    def custom_distribution_spherical(self, N, M_tot, inv_r_func, inv_theta_func, inv_phi_func):
        """
        Create initial positions and velocities of a custom distribution.

        Args:
            N (int): Number of particles.
            M_tot (float): Total mass.
            inv_r_func (function): Inverse function for radius.
            inv_theta_func (function): Inverse function for theta.
            inv_phi_func (function): Inverse function for phi.
        """
        self.N = N
        self.m = np.full(N, M_tot / N)
        
        # Sample random numbers for r, theta, and phi
        u_r = np.random.uniform(0, 1, N)
        u_theta = np.random.uniform(0, 1, N)
        u_phi = np.random.uniform(0, 1, N)
        
        # Generate positions using inverse functions
        self.r = inv_r_func(u_r)
        self.theta = inv_theta_func(u_theta)
        self.phi = inv_phi_func(u_phi)
        
        self.x = self.r * np.sin(self.theta) * np.cos(self.phi)
        self.y = self.r * np.sin(self.theta) * np.sin(self.phi)
        self.z = self.r * np.cos(self.theta)

    def custom_distribution_cyl(self, N, M_tot, inv_r_func, inv_theta_func, inv_z_func):
        """
        Create initial positions and velocities of a custom distribution.

        Args:
            N (int): Number of particles.
            M_tot (float): Total mass.
            inv_r_func (function): Inverse function for radius.
            inv_theta_func (function): Inverse function for theta.
            inv_phi_func (function): Inverse function for phi.
        """
        self.N = N
        self.m = np.full(N, M_tot / N)
        
        # Sample random numbers for r, theta, and phi
        u_r = np.random.uniform(0, 1, N)
        u_theta = np.random.uniform(0, 1, N)
        u_z = np.random.uniform(0, 1, N)
        
        # Generate positions using inverse functions
        self.r = inv_r_func(u_r)
        self.theta = inv_theta_func(u_theta)
        self.z = inv_z_func(u_z)
        
        
        self.x = self.r * np.cos(self.theta)
        self.y = self.r * np.sin(self.theta)
        self.z = self.z

    def custom_distribution_cartesian(self, N, M_tot, inv_x_func, inv_y_func, inv_z_func):
        """
        Create initial positions and velocities of a custom distribution.

        Args:
            N (int): Number of particles.
            M_tot (float): Total mass.
            inv_r_func (function): Inverse function for radius.
            inv_theta_func (function): Inverse function for theta.
            inv_phi_func (function): Inverse function for phi.
        """
        self.N = N
        self.m = np.full(N, M_tot / N)
        
        # Sample random numbers for r, theta, and phi
        u_x = np.random.uniform(0, 1, N)
        u_y = np.random.uniform(0, 1, N)
        u_z = np.random.uniform(0, 1, N)
        
        # Generate positions using inverse functions
        self.x = inv_x_func(u_x)
        self.y = inv_y_func(u_y)
        self.z = inv_z_func(u_z)


    def set_positions(self,N, M_tot, x, y, z):
        self.N = N
        self.m = np.full(N, M_tot / N)
        self.x = x
        self.y = y
        self.z = z



    # Include similar methods for other distribution models

    def initialize_velocities(self):
        self.vx = np.zeros(self.N)
        self.vy = np.zeros(self.N)
        self.vz = np.zeros(self.N)

    def set_velocities(self, vx, vy, vz):
        self.vx = vx
        self.vy = vy
        self.vz = vz


    def show_dist_sph(self, b=30):
        """
        Plot histograms of distributions.

        Args:
            b: number of bins
        """
        fig, axs = plt.subplots(1, 3, tight_layout=True, sharey=False)
        axs[0].set(title=fr'$r$ distribution', xlabel=r'$r$', ylabel=r'$p(r)$')
        axs[0].hist(self.r, bins=b, density=True, alpha=.5, color='royalblue', lw=0)
        axs[1].set(title=fr'$\theta$ distribution', xlabel=r'$\theta$', ylabel=r'$p(\theta)$')
        axs[1].hist(self.theta, bins=b, density=True, alpha=.5, color='royalblue', lw=0)
        axs[2].set(title=fr'$\phi$ distribution', xlabel=r'$\phi$', ylabel=r'$p(\phi)$')
        axs[2].hist(self.phi, bins=b, density=True, alpha=.5, color='royalblue', lw=0)

        return fig, axs

    def show_dist_cyl(self, b=30):
        """
        Plot histograms of distributions.

        Args:
            b: number of bins
        """
        fig, axs = plt.subplots(1, 3, tight_layout=True, sharey=False)
        axs[0].set(title=fr'$r$ distribution', xlabel=r'$r$', ylabel=r'$p(r)$')
        axs[0].hist(self.r, bins=b, density=True, alpha=.5, color='royalblue', lw=0)
        axs[1].set(title=fr'$\theta$ distribution', xlabel=r'$\theta$', ylabel=r'$p(\theta)$')
        axs[1].hist(self.theta, bins=b, density=True, alpha=.5, color='royalblue', lw=0)
        axs[2].set(title=fr'$z$ distribution', xlabel=r'$z$', ylabel=r'$p(z)$')
        axs[2].hist(self.z, bins=b, density=True, alpha=.5, color='royalblue', lw=0)

        return fig, axs

    def show_dist_cartesian(self, b=30):
        """
        Plot histograms of distributions.

        Args:
            b: number of bins
        """
        fig, axs = plt.subplots(1, 3, tight_layout=True, sharey=False)
        axs[0].set(title=fr'$r$ distribution', xlabel=r'$x$', ylabel=r'$p(x)$')
        axs[0].hist(self.x, bins=b, density=True, alpha=.5, color='royalblue', lw=0)
        axs[1].set(title=fr'$y$ distribution', xlabel=r'$y$', ylabel=r'$p(y)$')
        axs[1].hist(self.y, bins=b, density=True, alpha=.5, color='royalblue', lw=0)
        axs[2].set(title=fr'$z$ distribution', xlabel=r'$z$', ylabel=r'$p(z)$')
        axs[2].hist(self.z, bins=b, density=True, alpha=.5, color='royalblue', lw=0)
        
        return fig, axs

    def positions_2d(self, ax, axes='xy', l=10, a=1, s=1, bw=1, add_kde=False):
        """
        Plot 2D positions of particles.
        
        Args:
        - ax (matplotlib.axes._subplots.AxesSubplot): Axes object for plotting.
        - num (int): Time step index.
        - axes (str): Axes to plot ('xy', 'xz', or 'yz').
        - l (float): Limit of the domain.
        - a (float): Alpha value for scatter plot.
        """

        # Validate input axes
        valid_axes = {'x', 'y', 'z'}
        if not set(axes).issubset(valid_axes) or len(axes) != 2:
            print("Error: Axes must be specified as a combination of 'x', 'y', or 'z'.")
            return

        # Set plot title
        ax.set_title(f'{self.N} particles, mass = {round(np.sum(self.m), 2)}')

        # Map axis characters to corresponding data arrays
        axis_map = {'x': self.x, 'y': self.y, 'z': self.z}

        # Extract data for scatter plot
        x_data = axis_map[axes[0]]
        y_data = axis_map[axes[1]]
        
        if add_kde==True:
            xy = np.vstack([x_data, y_data])
            kde = gaussian_kde(xy, bw_method=bw)
            density = kde(xy)
        else:
            density='royalblue'
        
        
        ax.scatter(x_data, y_data, c=density, edgecolor='None', s=s, alpha=a)

        # Set x and y labels
        ax.set_xlabel(fr'${axes[0]}$')
        ax.set_ylabel(fr'${axes[1]}$')

        # Set limits, aspect ratio, and grid
        ax.set_xlim(-l, l)
        ax.set_ylim(-l, l)
        ax.set(aspect='equal')
        ax.grid()
        
        return ax
    

    def velocities_2d(self, ax, axes='xy', l=10, a=1, s=1, bw=1, add_kde=False):
        """
        Plot 2D positions of particles.
        
        Args:
        - ax (matplotlib.axes._subplots.AxesSubplot): Axes object for plotting.
        - num (int): Time step index.
        - axes (str): Axes to plot ('xy', 'xz', or 'yz').
        - l (float): Limit of the domain.
        - a (float): Alpha value for scatter plot.
        """

        # Validate input axes
        valid_axes = {'x', 'y', 'z'}
        if not set(axes).issubset(valid_axes) or len(axes) != 2:
            print("Error: Axes must be specified as a combination of 'x', 'y', or 'z'.")
            return

        # Set plot title
        ax.set_title(f'{self.N} particles, mass = {round(self.m[0]*self.N, 2)}')

        # Map axis characters to corresponding data arrays
        axis_map = {'x': self.vx, 'y': self.vy, 'z': self.vz}


        # Extract data for scatter plot
        x_data = axis_map[axes[0]]
        y_data = axis_map[axes[1]]

        if add_kde==True:
            xy = np.vstack([x_data, y_data])
            kde = gaussian_kde(xy, bw_method=bw)
            density = kde(xy)
        else:
            density='royalblue'
        
        
        ax.scatter(x_data, y_data, c=density, edgecolor='None', s=s, alpha=a)

        # Set x and y labels
        ax.set_xlabel(fr'${axes[0]}$')
        ax.set_ylabel(fr'${axes[1]}$')

        # Set limits, aspect ratio, and grid
        ax.set_xlim(-l, l)
        ax.set_ylim(-l, l)
        ax.set(aspect='equal')
        ax.grid()
        
        return ax


    def projections(self, l=10, a=1, s=1, bw=1, add_kde=False):
        fig, axs = plt.subplots(1,3, figsize = (11,3), tight_layout = True)
        self.positions_2d(axs[0], 'xy', l, a, s, bw, add_kde=add_kde)
        self.positions_2d(axs[1], 'xz', l, a, s, bw, add_kde=add_kde)
        self.positions_2d(axs[2], 'yz', l, a, s, bw, add_kde=add_kde)

        return fig, axs

    def v_projections(self, l=10, a=1, s=1, add_kde=False):
        fig, axs = plt.subplots(1,3, figsize = (11,11), tight_layout = True)
        self.velocities_2d(axs[0], 'xy', l, a, s, add_kde=add_kde)
        self.velocities_2d(axs[1], 'xz', l, a, s, add_kde=add_kde)
        self.velocities_2d(axs[2], 'yz', l, a, s, add_kde=add_kde)

        return fig, axs


    def positions_3d(self, ax, l=10, a=1, s=1):
        """
        Plot 3D positions of particles.

        Args:
            ax (mpl_toolkits.mplot3d.axes3d.Axes3D): 3D Axes object for plotting.
            num (int): Time step index.
            l (float): Limit of the domain.
            a (float): Alpha value for scatter plot.
        """
        ax.clear()

        # Set plot title and labels
        ax.set_title(f'{self.N} particles, mass = {round(self.m[0]*self.N, 2)}')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_zlabel(r'$z$')

        # Scatter plot of particle positions
        ax.scatter(self.x, self.y, self.z, color='royalblue', s=s, alpha=a)

        # Set limits and grid
        ax.set_xlim(-l, l)
        ax.set_ylim(-l, l)
        ax.set_zlim(-l, l)
        ax.set(aspect='equal')
        ax.grid()

        return ax


    def rho_mean_sphere(self, radius):
        """
        Calculates the mean density in a sphere of given radius
        """
        mask = self.r<radius
        M_r = np.sum(self.m[mask])
        V_r = 4/3*np.pi*radius**3
        rho = M_r/V_r

        return rho


    def write_sh1(self, output_file):
        """
        Write particle data to a file.

        Args:
            input_file (str): Path to the output file.
        """
        with open(output_file, 'w+') as f:
            f.write(str(self.N) + '\n')
            f.write('0' + '\n')
            for i in range(self.N):
                f.write(f'{self.m[i]} {self.x[i]} {self.y[i]} {self.z[i]} {self.vx[i]} {self.vy[i]} {self.vz[i]}\n')


    def write_tree(self, output_file):
        """
        Write particle data to a file.

        Args:
            input_file (str): Path to the output file.
        """
        with open(output_file, 'w+') as f:
            f.write(str(self.N) + '\n')
            f.write('3\n')  # Dimension
            f.write('0\n')  # Initial time
            for m_i in self.m:
                f.write(f'{m_i}\n')
            for i in range(self.N):
                f.write(f'{self.x[i]} {self.y[i]} {self.z[i]}\n')
            for j in range(self.N):
                f.write(f'{self.vx[j]} {self.vy[j]} {self.vz[j]}\n')



 


class plotter:
    """
    Plotting suite for SH1 Direct NBODY code
    by S.Mucedola
    """
    def __init__(self):
        """
        Initialize an instance of sh1_plotter.
        """
        # Initialize instance variables
        self.N = None  # Number of particles
        self.m = None  # Masses of particles
        self.t = None  # Time steps
        self.x, self.y, self.z = None, None, None  # Particle positions
        self.vx, self.vy, self.vz = None, None, None  # Particle velocities
        self.x_cm, self.y_cm , self.z_cm = None, None, None # com position
        self.r, self.v = None, None #distance from center of mass


    # -------------------- DATA ANALYSIS ----------------------- #

    def get_data_sh1(self, name_file):
        """
        Load data from a file and store it in the instance variables.
        
        Args:
        - name_file (str): Name of the file containing data.
        """
        # Read data from file
        f = np.array(open(name_file, 'r').readlines())
        
        # Extract information from data
        self.N = int(f[0])  # Number of particles
        s = self.N + 2  # Step size
        
        # Extract time steps
        self.t = f[1::s].astype(float)
        
        # Extract positions and velocities of particles
        mxv = np.array([[f[j+2::s][i].split(' ') for i in range(len(self.t))] for j in range(self.N)])
        mxv = mxv.astype(float)
        
        # Assign values to instance variables
        self.m = mxv[:, 0, 0]  # Masses
        x_ = mxv[:, :, 1:4]  # Positions
        v_ = mxv[:, :, 4:]  # Velocities
        self.x, self.y, self.z = [x_[:, :, i] for i in range(3)]
        self.vx, self.vy, self.vz = [v_[:, :, i] for i in range(3)]


        #calc CoM position and velocity
        self.x_cm, self.y_cm, self.z_cm = center_of_mass(self.x, self.y, self.z, self.m)
        self.vx_cm, self.vy_cm, self.vz_cm = center_of_mass(self.vx, self.vy, self.vz, self.m)
        
        
        #calc distance from com, velocity modulus
        self.r = distance(self.x, self.y, self.z, self.x_cm, self.y_cm, self.z_cm)
        self.v = distance(self.vx, self.vy, self.vz, self.vx_cm, self.vy_cm, self.vz_cm)


        # Print information about loaded data
        print('Loaded %.0e particles' % self.N)
        print('N particles, timesteps', self.x.shape)


    def get_data_tree(self, name_file):
        lines = pd.read_csv(name_file, names = list(range(3)), sep = '\s+')
        self.N = int(lines.iloc[0,0])
        s = self.N*3+3

        self.t = np.array([lines.iloc[2::s,0].copy().dropna()]).flatten()

        self.m = np.array([lines.iloc[i::s,0].copy().dropna() for i in range(3, self.N+3)])[:,0]
        x_ = np.array([lines.iloc[i::s].copy().dropna() for i in range(self.N+3, 2*self.N+3)])
        v_ = np.array([lines.iloc[i::s].copy().dropna() for i in range(2*self.N+3, 3*self.N+3)])
        
        self.x, self.y, self.z = [x_[:,:,i] for i in range(3)]
        self.vx, self.vy, self.vz = [v_[:,:,i] for i in range(3)]


        #calc CoM position and velocity
        M = np.sum(self.m)
        self.x_cm = np.sum(self.x * self.m[:, None], axis=0) / M
        self.y_cm = np.sum(self.y * self.m[:, None], axis=0) / M
        self.z_cm = np.sum(self.z * self.m[:, None], axis=0) / M
        self.vx_cm = np.sum(self.vx * self.m[:, None], axis=0) / M
        self.vy_cm = np.sum(self.vy * self.m[:, None], axis=0) / M
        self.vz_cm = np.sum(self.vz * self.m[:, None], axis=0) / M

        #calc distance from com, velocity modulus
        self.r = distance(self.x, self.y, self.z, self.x_cm, self.y_cm, self.z_cm)
        self.v = distance(self.vx, self.vy, self.vz, self.vx_cm, self.vy_cm, self.vz_cm)

        print('Loaded %.0e particles' %self.N)

    def convert_values(self, conversion_factors):
        for key in ['x', 'y', 'z', 'r']:
            setattr(self, key, getattr(self, key) * conversion_factors['pos'])
            
        for key in ['vx', 'vy', 'vz', 'v']:
            setattr(self, key, getattr(self, key) * conversion_factors['vel'])

        self.m *= conversion_factors['m']
        self.t *= conversion_factors['t']


    def get_log_tree(self, name_file):
        with open(name_file, "r") as file:
            lines = file.readlines()

        lines = [l for l in lines if l!='\n']
        lines = lines[3:]
        lines = [l.strip() for l in lines if not l.startswith('\t')]
        headers = lines[0].split()
        lines = [list(map(float,l.split())) for l in lines if not l.startswith('time')]

        d = pd.DataFrame(lines, columns=headers).to_dict()
        
        return d



    
    def particle_pos(self, particle):
        return self.x[particle,:], self.y[particle,:], self.z[particle,:] 
    def particle_vel(self, particle):
        return self.x[particle,:], self.y[particle,:], self.z[particle,:] 

    
    
    def rho_mean_sphere(self, radius):
        """
        Calculates the mean density in a sphere of given radius
        """
        rho = np.zeros(len(self.t))

        for num,_ in enumerate(self.t):
            mask = self.r[:,num]<radius
            M_r = np.sum(self.m[mask])
            V_r = 4/3*np.pi*radius**3
            rho[num] = M_r/V_r

        return rho
        
    def l_radii(self, qtile):
        """
        Returns the radius of a sphere with mass equal to a fraction of the total mass
        """
        L_radii = []
    
        # cycle over N_time_output
        for i in range(len(self.r[0])):
            radii_sorted = np.sort(self.r[:,i])
            mask = radii_sorted < np.quantile(radii_sorted, qtile)
            L_radii.append(np.max(radii_sorted[mask]))
            
        return np.array(L_radii)
        
        

    # ----------------- PLOTTING ROUTINES --------------------------- #


    def positions_2d(self, ax, num, axes='xy', l=10, b=0, a=1, s=10, lw=.5, bw=1, add_kde=False):
        """
        Plot 2D positions and trajectories of particles.
        
        Args:
        - ax (matplotlib.axes._subplots.AxesSubplot): Axes object for plotting.
        - num (int): Time step index.
        - axes (str): Axes to plot ('xy', 'xz', or 'yz').
        - l (float): Limit of the domain.
        - b (int): Length of trajectories.
        - a (float): Alpha value for scatter plot.
        """

        # Validate input axes
        valid_axes = {'x', 'y', 'z'}
        if not set(axes).issubset(valid_axes) or len(axes) != 2:
            print("Error: Axes must be specified as a combination of 'x', 'y', or 'z'.")
            return

        # Set plot title
        ax.set_title(f'{self.N} particles, mass = {round(np.sum(self.m), 2)}')

        # Map axis characters to corresponding data arrays
        axis_map = {'x': self.x, 'y': self.y, 'z': self.z}

        # Extract data for scatter plot
        x_data = axis_map[axes[0]][:,num]
        y_data = axis_map[axes[1]][:,num]
        
        if add_kde==True:
            xy = np.vstack([x_data, y_data])
            kde = gaussian_kde(xy, bw_method=bw)
            density = kde(xy)
        else:
            density='royalblue'
            
        
        # Set x and y labels
        ax.set_xlabel(fr'${axes[0]}$')
        ax.set_ylabel(fr'${axes[1]}$')

        # Plot scatter plot and trajectories
        ax.scatter(x_data, y_data, c=density, edgecolor='None', s=s, alpha=a)
        if b!=0:
            for i in range(self.N):
                b_ = round(b*num) 
                ax.plot(axis_map[axes[0]][i, num-b_:num], axis_map[axes[1]][i, num-b_:num], 
                        color='indianred', alpha=a, lw=lw)

        # Set limits, aspect ratio, and grid
        ax.set_xlim(-l, l)
        ax.set_ylim(-l, l)
        ax.set(aspect='equal')
        ax.grid()

        return ax
    
    
    def projections(self, num, l=10, b=0, a=1, s=1, lw=.5, bw=1, add_kde=False, axs=None):
        """
        Plot 2D projections of particles.

        Args:
            num (int): Time step index.
            l (float): Limit of the domain.
            b (int): Length of trajectories.
            a (float): Alpha value for scatter plot.
            s (float): Size of markers.
            lw (float): Width of trajectory lines.
            axs (list of matplotlib.axes._subplots.AxesSubplot, optional): List of Axes objects for plotting.
                If not provided, new figures will be created.

        Returns:
            None
        """
        if axs is None:
            fig, axs = plt.subplots(1, 3, figsize=(11, 11))
        else:
            fig = None

        self.positions_2d(axs[0], num, 'xy', l, b, a, s, lw, bw, add_kde=add_kde)
        self.positions_2d(axs[1], num, 'xz', l, b, a, s, lw, bw, add_kde=add_kde)
        self.positions_2d(axs[2], num, 'yz', l, b, a, s, lw, bw, add_kde=add_kde)

        return axs


    def positions_heatmap(self, num, axes='xy', x_min=-20, y_min=-20, l=10, bins=100, cmap='magma', ax=None):
        """
        Plot 2D heatmap of particle positions with custom range selection.

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot, optional): Axes object for plotting.
                If not provided, a new figure will be created.
            num (int): Time step index.
            axes (str): Axes to plot ('xy', 'xz', or 'yz').
            x_min (float): Minimum value of the x-axis range.
            y_min (float): Minimum value of the y-axis range.
            l (float): Length of the domain.
            bins (int or array_like or [int, int] or [array, array]): 
                Number of bins or bin edges.
            cmap (str or Colormap): Colormap for the heatmap.

        Returns:
            matplotlib.axes._subplots.AxesSubplot: Axes object with heatmap plot.
        """

        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 5))
        else:
            fig = None

        # Validate input axes
        valid_axes = {'x', 'y', 'z'}
        if not set(axes).issubset(valid_axes) or len(axes) != 2:
            print("Error: Axes must be specified as a combination of 'x', 'y', or 'z'.")
            return

        # Set plot title
        ax.set_title(f'{self.N:.0e} particles, mass = {self.m[0]*self.N:.1e}, time = {self.t[num]:.1f}')

        # Map axis characters to corresponding data arrays
        axis_map = {'x': self.x, 'y': self.y, 'z': self.z}

        # Extract data for heatmap
        x_data = axis_map[axes[0]][:, num]
        y_data = axis_map[axes[1]][:, num]

        range_x = [x_min, x_min + l]
        range_y = [y_min, y_min + l]

        # Plot heatmap with custom range selection
        ax.hist2d(x_data, y_data, bins=bins, range=[range_x, range_y], cmap=cmap)

        # Set x and y labels
        ax.set_xlabel(fr'${axes[0]}$')
        ax.set_ylabel(fr'${axes[1]}$')

        # Set aspect ratio
        ax.set_aspect('equal')

        return ax


    def projections_heatmap(self, num, x_min=-10, y_min=-10, z_min=-10, l=10, bins=100, cmap='magma', axs=None):
        """
        Plot 2D projections of particles using heatmaps.

        Args:
            num (int): Time step index.
            x_min (float): Minimum value of the x-axis range.
            y_min (float): Minimum value of the y-axis range.
            z_min (float): Minimum value of the z-axis range.
            l (float): Length of the domain.
            bins (int or array_like): Number of bins or bin edges for the histograms.
            cmap (str or Colormap): Colormap for the heatmaps.
            axs (list of matplotlib.axes._subplots.AxesSubplot, optional): List of Axes objects for plotting.
                If not provided, new figures will be created.

        Returns:
            list of matplotlib.axes._subplots.AxesSubplot: List of Axes objects with heatmap plots.
        """
        if axs is None:
            fig, axs = plt.subplots(1, 3, figsize=(11, 11), tight_layout=True)
        else:
            fig = None


        # Plot each projection heatmap
        self.positions_heatmap(num, 'xy', x_min, y_min, l, bins, cmap, ax=axs[0])
        self.positions_heatmap(num, 'xz', x_min, z_min, l, bins, cmap, ax=axs[1])
        self.positions_heatmap(num, 'yz', y_min, z_min, l, bins, cmap, ax=axs[2])


        return fig, axs



    def positions_3d(self, ax, num, l=10, b=0, a=1, s=10, lw=.5):
        """
        Plot 3D positions and trajectories of particles.

        Args:
        - ax (mpl_toolkits.mplot3d.axes3d.Axes3D): 3D Axes object for plotting.
        - num (int): Time step index.
        - l (float): Limit of the domain.
        - b (int): Length of trajectories.
        - a (float): Alpha value for scatter plot.
        """
        ax.clear()

        # Set plot title and labels
        ax.set_title(f'{self.N:n} particles, mass = {self.m[0]*self.N:n}, time = {self.t[num]:.1f}')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_zlabel(r'$z$')

        # Scatter plot of particle positions
        ax.scatter(self.x[:, num], self.y[:, num], self.z[:, num], color='royalblue', s=s, alpha=a)

        # Plot trajectories
        if b!=0:
            for i in range(self.N):
                b_ = round(b*num)
                ax.plot(self.x[i, num-b_:num], self.y[i, num-b_:num], self.z[i, num-b_:num],
                        color='indianred', alpha=a, lw=lw)

        # Set limits and grid
        ax.set_xlim(-l, l)
        ax.set_ylim(-l, l)
        ax.set_zlim(-l, l)
        ax.set(aspect='equal')
        ax.grid()

        return ax
    
    def heatmap_3d(self, num=0, s=1, alpha=.5, xl=10, yl=10, zl=5):
        fig = go.Figure()

        x = self.x[:,num]
        y = self.y[:,num]
        z = self.z[:,num]


        # Add scatter plot trace with density as color variable
        scatter = go.Scatter3d(
        x=x,y=y,z=z,
        mode='markers',
        marker=dict(
            size=s,
            opacity=alpha
        )
        )


        # Add trace to the figure
        fig.add_trace(scatter)

        # Update layout options for better visualization (optional)
        fig.update_scenes(
        dict(
            xaxis=dict(range=[-xl, xl]),   # Set xlim
            yaxis=dict(range=[-yl, yl]),   # Set ylim
            zaxis=dict(range=[-zl, zl])    # Set zlim
        )
        )

        # Show the figure
        fig.show()
    

    def plot_mass_distribution(self, num, x=None, y=None, z=None, labels=None, bins=500, xlim=(-20, 20)):
        """
        Plot mass distribution along specified axes with medians.

        Args:
            num (int): Time step index.
            x, y, z (array-like or None, optional): Data arrays for the x, y, and z coordinates.
                If None, uses self.x, self.y, and self.z as default.
            labels (list of str, optional): List of labels for the data arrays.
                If None, uses ['x', 'y', 'z'] as default.
            bins (int): Number of bins for the histograms.
            xlim (tuple): Limits for the x-axis.

        Returns:
            axis
        """
        fig, ax = plt.subplots()

        if x is None and y is None and z is None:
            x = self.x
            y = self.y
            z = self.z
        if labels is None:
            labels = ['x', 'y', 'z']



        counts_x, bins_x, _ = ax.hist(x[:, num], bins=bins, histtype='bar', edgecolor='g', color='g', alpha=.2,
                                    label=r'along $x$')
        counts_y, bins_y, _ = ax.hist(y[:, num], bins=bins, histtype='bar', edgecolor='orange', color='orange',
                                    alpha=.2, label=r'along $y$')
        counts_z, bins_z, _ = ax.hist(z[:, num], bins=bins, histtype='bar', edgecolor='royalblue', color='royalblue',
                                    alpha=.2, label=r'along $z$')
        ax.set_title('Mass distribution origin reference frame')
        ax.set_xlabel('Distance from origin (pc)')
        ax.set_ylabel(r'Mass ($\mathrm{M_\odot}$)')
        ax.set_xlim(*xlim)

        # Calculate medians
        median_x = np.median(x[:, num])
        median_y = np.median(y[:, num])
        median_z = np.median(z[:, num])

        # Plot medians
        ax.axvline(median_x, ymin=0, ymax=np.max(counts_x) / ax.get_ylim()[1], color='g', lw=.8,
                label=fr'median $x = {median_x:.1f}$ pc')
        ax.axvline(median_y, ymin=0, ymax=np.max(counts_y) / ax.get_ylim()[1], color='orange', lw=1,
                label=fr'median $y = {median_y:.1f}$ pc')
        ax.axvline(median_z, ymin=0, ymax=np.max(counts_z) / ax.get_ylim()[1], color='royalblue', lw=1,
                label=fr'median $z = {median_z:.1f}$ pc')

        ax.legend()
        
        return ax


    def preview_animation(self, step=5, l=100, lz=10):
        """
        Create a 3D scatter plot animation of particle positions.

        Args:
            step (int): Time step interval for animation.
            l (float): Range limit for x and y axes.
            lz (float): Range limit for z axis.

        Returns:
            plotly.graph_objects.Figure: 3D scatter plot animation figure.
        """
        # Slicing
        x = self.x[:, ::step]
        y = self.y[:, ::step]
        z = self.z[:, ::step]

        n_particles = x.shape[0]
        times = x.shape[1]

        # Create a long-format DataFrame
        data = {
            'x': x.flatten(),
            'y': y.flatten(),
            'z': z.flatten(),
            'particle': np.repeat(np.arange(n_particles), times),
            'time': np.tile(np.arange(times), n_particles)
        }

        df = pd.DataFrame(data)

        # Create the 3D scatter plot animation
        fig = px.scatter_3d(df, 
                            x='x', 
                            y='y', 
                            z='z',
                            animation_frame='time',
                            range_x=[-l, l], 
                            range_y=[-l, l], 
                            range_z=[-lz, lz])

        fig.update_traces(marker=dict(size=1, opacity=0.2))  # Adjust size and opacity
        return fig



    ### animations
    def animate_2d(self, fig, ax, axes='xy', f_t=1, l=10, b=0, a=1, s=10, lw=.5, frame_skip=1, interval=1):
        """
        Animate 2D positions and trajectories of particles over time.
        
        Args:
        - ax (matplotlib.axes._subplots.AxesSubplot): Axes object for plotting.
        - axes (str): Axes to plot ('xy', 'xz', or 'yz').
        - f_t (float): fraction of total time.
        - l (float): Limit of the domain.
        - b (int): Length of trajectories.
        
        Returns:
        - animation.FuncAnimation: Animation object.
        """
        # Define animation function
        def update(num):
            ax.clear()
            self.positions_2d(ax, num=num, axes=axes, l=l, b=b, a=a, s=s, lw=lw)
        
        # Create animation
        num_frames = round(f_t*len(self.t))
        ani_2d = animation.FuncAnimation(fig, func=update, frames=tqdm(range(num_frames)[::frame_skip]), interval=interval)
        return ani_2d
    
    
    def animate_projections(self, fig, axs, f_t=1, l=10, b=0, a=1, s=10, lw=.5, frame_skip=1, interval=1):
        """
        Animate projections of particles over time.

        Args:
        - fig (matplotlib.figure.Figure): Figure object for the entire plot.
        - axs (list of matplotlib.axes._subplots.AxesSubplot): List of Axes objects for plotting projections.
        - f_t (float): Fraction of total time.
        - l (float): Limit of the domain.
        - b (int): Length of trajectories.
        - a (float): Alpha value for scatter plot.
        - s (float): Size of markers.
        - lw (float): Width of trajectory lines.
        - frame_skip (int): Number of frames to skip.
        - interval (int): Interval between frames in milliseconds.

        Returns:
        - animation.FuncAnimation: Animation object.
        """
        # Define animation function
        def update(num):
            for ax, axes in zip(axs, ['xy', 'xz', 'yz']):
                ax.clear()
                self.positions_2d(ax, num=num, axes=axes, l=l, b=b, a=a, s=s, lw=lw)

        # Create animation
        num_frames = round(f_t * len(self.t))
        ani_projections = animation.FuncAnimation(fig, func=update, frames=tqdm(range(num_frames)[::frame_skip]),
                                                  interval=interval)
        return ani_projections



    def animate_3d(self, fig, ax, f_t=1, l=10, b=0, a=1, s=10, lw=.5, frame_skip=1, interval=1):
        """
        Animate 3D positions and trajectories of particles over time.
        
        Args:
        - ax (mpl_toolkits.mplot3d.axes3d.Axes3D): 3D Axes object for plotting.
        - f_t (float): fraction of total time.
        - l (float): Limit of the domain.
        - b (int): Length of trajectories.
        
        Returns:
        - animation.FuncAnimation: Animation object.
        """
        # Define animation function
        def update(num):
            ax.clear()
            self.positions_3d(ax, num=num, l=l, b=b, a=a, s=s, lw=lw)
            scatter = ax.scatter(self.x[:, num], self.y[:, num], self.z[:, num], c='royalblue', s=s, alpha=a)
            return scatter

        # Create animation
        num_frames = round(f_t*len(self.t))
        ani_3d = animation.FuncAnimation(fig, func=update, frames=tqdm(range(num_frames)[::frame_skip]), interval=interval, blit=True)
        return ani_3d
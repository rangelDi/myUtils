import numpy as np
import matplotlib.pyplot as plt
import sys

""" MISC """
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def exclusiveExistence(a: float, b: float) -> bool:
    return (np.isnan(a) and not np.isnan(b)) or (np.isnan(b) and not np.isnan(a))

""" MAIN CLASS """
class AdvectionDiffusion2D:

    def __init__(self, x_range: list[float], N_x: float, y_range: list[float], N_y: float,
                 gamma: float, rho: float, u: float, v: float, boundaryFlux: list[float],
                 boundaryPhi: list[float]):

        """ DOMAIN SETUP """
        self.x_min: float = x_range[0]
        self.x_max: float = x_range[1]

        self.y_min: float = y_range[0]
        self.y_max: float = y_range[1] 

        self.N_x: float = N_x
        self.N_y: float = N_y

        """ MESH SETUP """
        self.L_y: float = self.y_max - self.y_min 
        self.L_x: float = self.x_max - self.x_min 
        self.dx: float  = self.L_x/self.N_x 
        self.dy: float  = self.L_x/self.N_y 
        print(f"dy: {self.dy}")
        print(f"dx: {self.dx}")

        self.x_gridlines: np.ndarray = np.linspace(self.x_min, self.x_max, self.N_x + 1)
        self.y_gridlines: np.ndarray = np.linspace(self.y_min, self.y_max, self.N_y + 1) 

        self.x_c: np.ndarray = np.zeros((self.N_x + 2, self.N_y + 2))
        self.y_c: np.ndarray = np.zeros((self.N_x + 2, self.N_y + 2))

        for i in range(1, self.N_x+1):
            self.x_c[i] = self.x_gridlines[i-1] + 0.5*self.dx
        for i in range(1, self.N_y+1):
            self.y_c[i] = self.y_gridlines[i-1] + 0.5*self.dy

        self.x_c[0] = self.x_gridlines[0]
        self.y_c[0] = self.y_gridlines[0]
        self.x_c[-1] = self.x_gridlines[-1]
        self.y_c[-1] = self.y_gridlines[-1]

        self.x_c = self.x_c.T
        self.phi = np.zeros((self.N_x, self.N_y))

        """ PROPERY SETUP """
        self.gamma: float = gamma
        self.rho: float = rho if not np.isnan(rho) else 0
        self.u: float = u
        self.v: float = v

        self.F_x: float = self.rho * self.u * self.dy
        self.F_y: float = self.rho * self.v * self.dx
        self.D_x: float = self.gamma * self.dy/self.dx
        self.D_y: float = self.gamma * self.dx/self.dy
        self.Pe_x: float = self.F_x/self.D_x
        self.Pe_y: float = self.F_y/self.D_y 
        print(f"Pe_x: {self.Pe_x}")
        print(f"Pe_y: {self.Pe_y}")

        """ BOUNDARY CONDITION SETUP """
        numberOfBoundaries: int = 4
        if len(boundaryFlux) != numberOfBoundaries or len(boundaryPhi) != numberOfBoundaries:
            raise ValueError(f"Both 'boundaryFlux' and 'boundaryPhi' must have length {numberOfBoundaries}")
        
        for i in range(0,numberOfBoundaries):
            if not exclusiveExistence(boundaryFlux[i], boundaryPhi[i]):
                raise ValueError(f"Every boundary may have only a prescribed phi or a prescribed flux.")

        self.flux_n: float = boundaryFlux[0]
        self.flux_e: float = boundaryFlux[1]
        self.flux_s: float = boundaryFlux[2]
        self.flux_w: float = boundaryFlux[3]

        self.phi_n: float = boundaryPhi[0]
        self.phi_e: float = boundaryPhi[1]
        self.phi_s: float = boundaryPhi[2]
        self.phi_w: float = boundaryPhi[3]
    
        advectiveConditions: bool = (self.u != 0.0) or (self.v != 0.0)
        self.prescribedFlux: bool = not(np.isnan(self.flux_n) and np.isnan(self.flux_e)
                                   and np.isnan(self.flux_s) and np.isnan(self.flux_w))
        
        if advectiveConditions and self.prescribedFlux:
            raise ValueError("Prescribed flux not implemented for advective conditions")

        if advectiveConditions and np.isnan(self.rho):
            raise ValueError("Rho must be provided for advective conditions")

        """ SET INTERNAL COEFFICIENTS """
        self.A_w: np.ndarray = np.ones((self.N_x, self.N_y)) * max([ self.F_x, self.D_x + self.F_x/2, 0]) 
        self.A_e: np.ndarray = np.ones((self.N_x, self.N_y)) * max([-self.F_x, self.D_x - self.F_x/2, 0])
        self.A_s: np.ndarray = np.ones((self.N_x, self.N_y)) * max([ self.F_y, self.D_y + self.F_y/2, 0])
        self.A_n: np.ndarray = np.ones((self.N_x, self.N_y)) * max([-self.F_y, self.D_y - self.F_y/2, 0])
        self.q_p: np.ndarray = np.zeros((self.N_x, self.N_y))
        self.q_u: np.ndarray = np.zeros((self.N_x, self.N_y))
        self.A_p: np.ndarray = 0

        """ SET BOUNDARY CONDITIONS """
        if self.prescribedFlux:
            self.setDiffusionBCs()
        else:
            self.setAdvectionDiffusionBCs()
        
        """ PREPARE SYSTEM """
        self.A_p = self.A_w + self.A_e + self.A_s + self.A_n - self.q_p;

    """ BOUNDARY CONDITIONS """
    def setAdvectionDiffusionBCs(self):
        if(self.Pe_x > 2):
            self.A_w[0,:]   = 0
            self.q_p[0,:]  += -(2*self.D_x+self.F_x)
            self.q_u[0,:]  += (2*self.D_x+self.F_x)*self.phi_w

            self.A_e[self.N_x-1,:]  = 0  
            self.q_p[self.N_x-1,:] += -(2*self.D_x)  
            self.q_u[self.N_x-1,:] += (2*self.D_x)*self.phi_e  
        elif(self.Pe_x < -2):
            self.A_w[0,:]   = 0
            self.q_p[0,:]  += -(2*self.D_x)
            self.q_u[0,:]  += (2*self.D_x)*self.phi_w

            self.A_e[self.N_x-1,:]  = 0  
            self.q_p[self.N_x-1,:] += self.F_x-2*self.D_x  
            self.q_u[self.N_x-1,:] += -(self.F_x-2*self.D_x)*self.phi_e  
        else:
            self.A_w[0,:]   = 0
            self.q_p[0,:]  += -(2*self.D_x+self.F_x)
            self.q_u[0,:]  += (2*self.D_x+self.F_x)*self.phi_w

            self.A_e[self.N_x-1,:]  = 0  
            self.q_p[self.N_x-1,:] += self.F_x-2*self.D_x  
            self.q_u[self.N_x-1,:] += -(self.F_x-2*self.D_x)*self.phi_e  

        if(self.Pe_y > 2):
            self.A_s[:,0]   = 0
            self.q_p[:,0]  += -(2*self.D_y+self.F_y)
            self.q_u[:,0]  += (2*self.D_y+self.F_y)*self.phi_s

            self.A_n[:,self.N_y-1]  = 0  
            self.q_p[:,self.N_y-1] += -(2*self.D_y)  
            self.q_u[:,self.N_y-1] += (2*self.D_y)*self.phi_n  
        elif(self.Pe_y < -2):
            self.A_s[:,0]   = 0
            self.q_p[:,0]  += -(2*self.D_y)
            self.q_u[:,0]  += (2*self.D_y)*self.phi_s

            self.A_n[:,self.N_y-1]  = 0  
            self.q_p[:,self.N_y-1] += self.F_y-2*self.D_y  
            self.q_u[:,self.N_y-1] += -(self.F_y-2*self.D_y)*self.phi_n  
        else:
            self.A_s[:,0]   = 0
            self.q_p[:,0]  += -(2*self.D_y+self.F_y)
            self.q_u[:,0]  += (2*self.D_y+self.F_y)*self.phi_s

            self.A_n[:,self.N_y-1]  = 0  
            self.q_p[:,self.N_y-1] += self.F_y-2*self.D_y  
            self.q_u[:,self.N_y-1] += -(self.F_y-2*self.D_y)*self.phi_n  

    """ SET DIFFUSION ONLY BOUNDARY CONDITIONS """
    def setDiffusionBCs(self):
        if(np.isnan(self.flux_n)):
            self.A_n[:,self.N_y-1]  = 0  
            self.q_p[:,self.N_y-1] += self.F_y-2*self.D_y  
            self.q_u[:,self.N_y-1] += -(self.F_y-2*self.D_y)*self.phi_n  
        else: # flux_n is applied
            self.A_n[:,self.N_y-1]  = 0  
            self.q_p[:,self.N_y-1] += self.F_y  
            self.q_u[:,self.N_y-1] += self.flux_n  

        if(np.isnan(self.flux_e)):
            self.A_e[self.N_x-1,:]  = 0  
            self.q_p[self.N_x-1,:] += self.F_x-2*self.D_x  
            self.q_u[self.N_x-1,:] += -(self.F_x-2*self.D_x)*self.phi_e  
        else: # flux_e is applied
            self.A_e[self.N_x-1,:]  = 0  
            self.q_p[self.N_x-1,:] += self.F_x  
            self.q_u[self.N_x-1,:] += self.flux_e  

        if(np.isnan(self.flux_s)):
            self.A_s[:,0]   = 0
            self.q_p[:,0]  += -(2*self.D_y+self.F_y)
            self.q_u[:,0]  += (2*self.D_y+self.F_y)*self.phi_s
        else: # flux_s is applied
            self.A_s[:,0]   = 0
            self.q_p[:,0]  += -self.F_y
            self.q_u[:,0]  += self.flux_s

        if(np.isnan(self.flux_w)):
            self.A_w[0,:]   = 0
            self.q_p[0,:]  += -(2*self.D_x+self.F_x)
            self.q_u[0,:]  += (2*self.D_x+self.F_x)*self.phi_w
        else: # flux_w is applied
            self.A_w[0,:]   = 0
            self.q_p[0,:]  += -self.F_x
            self.q_u[0,:]  += self.flux_w

    def TDMA(self, leftDiagonal_in, centralDiagonal_in, rightDiagonal_in, righthandSide_in):
        leftDiagonal = leftDiagonal_in.copy()
        centralDiagonal = centralDiagonal_in.copy()
        rightDiagonal = rightDiagonal_in.copy()
        righthandSide = righthandSide_in.copy()

        N: int = len(centralDiagonal)  # size of system
        x: np.ndarray = np.zeros(N)  # solution column vector

        # leftDiagonal is a vector of length N but entry 1 is 0
        # rightDiagonal is a vector of length N but entry N is 0

        # FORWARD PASS
        for i in range(1, N):  # [1, N-1] inclusive
            l: float = leftDiagonal[i] / centralDiagonal[i-1]  # elimination coefficient
            centralDiagonal[i] = centralDiagonal[i] - l * rightDiagonal[i-1]
            righthandSide[i] = righthandSide[i] - l * righthandSide[i-1]

        # FIRST SOLUTION
        # centralDiagonal[N-1] * x[N-1] = righthandSide[N-1]
        x[N-1] = righthandSide[N-1] / centralDiagonal[N-1]

        # BACKWARD PASS
        # centralDiagonal[i] * x[i] + rightDiagonal[i] * x[i+1] = righthandSide[i]
        for i in range(N-2, -1, -1):  # [N-2, 0] inclusive, backwards step
            x[i] = (righthandSide[i] - rightDiagonal[i] * x[i+1]) / centralDiagonal[i]

        return x

    def iterative2DSolver(self, A_n, A_s, A_w, A_e, A_p, q_u, N_x, N_y, iterations, tolerance=1e-4):
        Q = np.zeros((N_x, N_y))
        phi = np.zeros((N_x, N_y))
        residual = 0
        normalizedResidual = 0

        for iteration in range(iterations):
            for j in range(1, N_y-1):
                Q[:,j] = q_u[:,j] + A_n[:,j] * phi[:,j+1] + A_s[:,j] * phi[:,j-1]

            Q[:,0] = q_u[:,0] + A_n[:,0] * phi[:,1] + A_s[:,0] * phi[:,0]  # <- set phi(:,0) = phi(:,0)
            Q[:,N_y-1] = q_u[:,N_y-1] + A_n[:,N_y-1] * phi[:,N_y-1] + A_s[:,N_y-1] * phi[:,N_y-2]  # <- set phi(:,N_y) = phi(:, N_y-1)

            for j in range(N_y):
                phi[:,j] = self.TDMA(-A_w[:,j], A_p[:,j], -A_e[:,j], Q[:,j])  # <- solve each row as if it was its own 1D problem

            residual = 0  # <- residual should ideally be 0
            for i in range(1, N_x-1):
                for j in range(1, N_y-1):
                    residual = residual + np.abs(-A_e[i,j] * phi[i+1,j] - A_w[i,j] * phi[i-1,j] - A_n[i,j] * \
                        phi[i,j+1] - A_s[i,j] * phi[i,j-1] + A_p[i,j] * phi[i,j] - q_u[i,j])

            if iteration == 0:
                normalizedResidual = residual

            divByZeroGuard = 1e-20
            RSM = residual / (normalizedResidual + divByZeroGuard)

            if RSM < tolerance:
                break

        return phi

    def solve(self, iterations = 1000):

        self.phi = self.iterative2DSolver(self.A_n, self.A_s, self.A_w, self.A_e, self.A_p, self.q_u, self.N_x, self.N_y, iterations);

        # include boundary phi
        phi_all = np.zeros((self.N_x+2, self.N_y+2));
        phi_all[1:self.N_x+1, 1:self.N_y+1] = self.phi;
        self.phi = phi_all

        # northern boundary
        if np.isnan(self.flux_n):
            self.phi[:, -1] = self.phi_n
        else:
            self.phi[:, -1] = self.phi[:, -2] - self.flux_n * self.dy / (2 * self.gamma)

        # eastern boundary
        if np.isnan(self.flux_e):
            self.phi[-1, :] = self.phi_e
        else:
            self.phi[-1, :] = self.phi[-2, :] - self.flux_e * self.dx / (2 * self.gamma)

        # southern boundary
        if np.isnan(self.flux_s):
            self.phi[:, 0] = self.phi_s
        else:
            self.phi[:, 0] = self.phi[:, 1] + self.flux_s * self.dy / (2 * self.gamma)

        # western boundary
        if np.isnan(self.flux_w):
            self.phi[0, :] = self.phi_w
        else:
            self.phi[0, :] = self.phi[1, :] + self.flux_w * self.dx / (2 * self.gamma)

        self.phi = self.phi.T
        return self.x_c, self.y_c, self.phi

    def plot(self, xlabel = 'x', ylabel = 'y', title = r'$\phi$ Distribution', cbarLabel = r"Units of $\phi$"):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(self.x_c, self.y_c, self.phi, cmap='viridis', rasterized = False)
        
        ax.set_proj_type('ortho')
        ax.view_init(90, 270)
        cbar = fig.colorbar(surf, ax=ax)
        cbar.set_label(cbarLabel)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        ax.set_zticks([])
        plt.title(title)


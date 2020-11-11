#-----------------------------------------------------------------------------------
# Reference:  P. Raanes
# https://github.com/nansencenter/DAPPER/blob/master/da_methods/da_methods.py
# Reference: Ceci Dip LIAMG
# https://github.com/groupeLIAMG/EnKF/blob/master/ENKF.py
# Reference: Evensen, Geir. (2009):
# "The ensemble Kalman filter for combined state and parameter estimation."
#-----------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
# Ensemble Kalman Filter class by Thomas Beraud, 2019
# INRS ETE Quebec
#-----------------------------------------------------------------------------------

#------- Import -------#

from random import gauss
import numpy as np
import scipy as sp
from scipy.linalg import solve
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy.random import multivariate_normal
from scipy.sparse import csr_matrix

""" --------------   Importing Ensemble dataset   --------------

The ensemble matrix should be formated as each row
contained a parameter and each column a simulation

If you run 100 simulations with 10 000 grid parameters
to estimate, the ensemble matrix is of shape :

(10 000, 100)
"""

""" --------------   Importing Observations Dataset   --------------

The ensemble matrix should be formated as each row
contained one observation and each column a simulation

To use an Ensemble Smoother you should update the model by
assimilating all the data in one run. So you have to stack
all the observed data in one matrix

If you run 100 simulations, with 20 observations parameters
at each timestep and 50 timestep, the observation matrix is of shape :

1000 observations = 20 * 50

(1000, 100)

Observation matrix :

|   observations at 1 time-step  |
|   observations at 2 time-ste   |
|               .                |
|               .                |
| observations at last time-step |
"""


""" --------------   Importing Observations Localisation   --------------


Each observation have to be localized in the grid parameters we are
evaluating. For each observation we have to had the number of the
corresponding element in the grid.

The localisation vector have to be the same length as the observation matrix.
This is possible to look at different point at each time-step. To allow this
flexibilty user have to defined the complete vector himself. If this is only
the same point that are sampled in the gridat each time-step, the user can
simply repeat the observations vector for the number of timestep.

In background, a sparse observation matrix will be build to allow computation
between the parameters matrix and only the sampled points.

If you run 100 simulations, with 20 observations parameters
at each timestep and 50 timestep, the localisation vector is of shape :

(1000, 1)

Observation matrix :

|   localisation of observations at 1 time-step  |
|   localisation of observations at 2 time-ste   |
|               .                |
|               .                |
| localisation of observations at last time-step |
"""


def Measurements(  measurements_vector,
                   measurements_error_percent,
                   number_of_ensemble,
                   alpha=1,
                   show_stats=False,
                   show_cov_matrix=False,
                  ):

        d = measurements_vector
        N = number_of_ensemble
        m = measurements_vector.shape[0]
        D = np.zeros((N,m))
        Cov = np.zeros((m,m))
        var_temp = np.zeros((m,N))

        #D = multivariate_normal([0]*test.p_, test.Cd_, test.N_,check_valid='warn')

        np.random.seed( 42 )

        for simu in range(N):
            for elt in range(len(d)):
                temp = gauss(d[elt],(measurements_error_percent)*d[elt])
                if temp >= 0 :
                    D[simu][elt] = temp
                else :
                    D[simu][elt] = 0


        D = np.matrix(D.T)
        #mean_D = np.mean(D,1)
        mean_D = np.matrix(d).T

        anomalies = D - mean_D
        Cov = anomalies @ anomalies.T / (N - 1)


        if show_stats:
            print("#----- Measurements statistics -----#\n"+
            "Number of measures : %.0f" % d.shape[0] +"\n"+
            "Mean of measures : %.6f" % np.mean(d) +"\n"+
            "Variance of measures : %.6f" % np.var(d,dtype=np.float64) +"\n"+
            "Standard Deviation of measures : %.4f" % np.sqrt(np.var(d,dtype=np.float64)) +"\n"+
            "#----------------------------------#\n")

        if show_cov_matrix:

            plt.figure(figsize=(12,8))
            plt.imshow(D,aspect='auto')
            plt.colorbar()
            plt.title('Data Measurements Matrix after perturbation')
            plt.show()


            plt.figure(figsize=(12,8))
            plt.imshow(Cov)
            plt.title('Covariance Measurements Matrix')
            plt.colorbar()
            plt.show()



        return (Cov,D)

def Observation_Matrix(
               vector_of_cell_observation,
               number_of_total_grid_cell,
               show_observation_matrix=False,
            ):

        """
        #Test   :
        H = Observation_Matrix(np.array([12,36,23,45,67,78,87]),100,True)
        """

        number_of_observations = vector_of_cell_observation.shape[0]

        # Using sparse matrix to allow efficient computation,
        # almost the whole matrix is nil

        # Row of observation
        row = vector_of_cell_observation

        # Column of observation
        col = np.zeros(number_of_observations)
        for elt in range(number_of_observations):
            col[elt] = elt

        # Observed point are simply equal to one in matrix product
        observed_point = np.ones(number_of_observations)

        H = csr_matrix((observed_point, (row, col)), shape=(number_of_total_grid_cell, number_of_observations))

        return (H)

def cell_index_to_position(cell_number,nx=41,ny=41,nz=51):
    return (cell_number%nx),(int(cell_number/nx))

def cell_position_to_index(x,y,z,nx=41,ny=41,nz=51):
    return (z*(nx*ny)+y*nx+x)

def pseudo_inverse(Cd,Cdd,alpha=1,method='subspace'):
    """ --------------   From the paper of Emerick, 2012   --------------
                        DOI: 10.1007/s10596-012-9275-5

    The two covariances matrices are given as parameters
    Cd is the covariance matrix of observed data measurements
    Cdd is the auto-covariance matrix of predicted data

    Returns the pseudo-inverse of the covaraiance matrix C
    """
    if method=='tsvd':







        C = alpha*Cd + Cdd

        U,S,V = np.linalg.svd(C)

        pseudo = V.T@np.linalg.inv(np.diag(S))@U.T
        """

        alpha_Cd = alpha*Cd

        regul_term = 0.5 #1.0e-5
        while(np.linalg.det(alpha_Cd)==0):
            print("Regul term : ",regul_term)
            alpha_Cd = alpha_Cd + np.eye(alpha_Cd.shape[0])*regul_term
            regul_term += 0.1




        lower_Cd = np.linalg.cholesky(alpha_Cd)
        upper_Cd = lower_Cd.T
        inv_lower_Cd = np.linalg.solve(lower_Cd,np.eye(lower_Cd.shape[0]))
        inv_upper_Cd = np.linalg.solve(upper_Cd,np.eye(upper_Cd.shape[0]))

        C = inv_lower_Cd@Cdd@inv_upper_Cd + np.eye(Cd.shape[0])

        U,S,V = np.linalg.svd(C,full_matrices=False)

        pseudo = inv_upper_Cd@V.T@np.linalg.solve(np.diag(S),np.eye(S.shape[0]))@U.T@inv_lower_Cd
        """
    if method=='subspace':
        print('not implemented')





    else :
        "Wrong method name, you should choose 'tsvd' (if m<N) or 'subspace'"
    return pseudo


class EnsembleKalmanFilter_Emerick():


    def __init__(self,
                   ensemble_matrix,
                   measurements_vector,
                   observations_vector,
                   measurements_error_percent,
                   alpha,
                   show_parameters=False,
                  ):


        self.m_ = ensemble_matrix.shape[0]
        self.N_ = ensemble_matrix.shape[1]
        self.p_ = measurements_vector.shape[0]
        self.Ef_ = ensemble_matrix
        self.Ea_ = np.zeros((ensemble_matrix.shape))
        self.y_ = measurements_vector
        self.hE_ = np.zeros((self.p_,self.N_))
        self.alpha_ = alpha
        self.Cd_ = np.zeros((self.m_,self.m_))
        self.measurements_error_percent_ = measurements_error_percent
        self.display_ = show_parameters

        if show_parameters :
            print("#----- Ensemble Kalman Filter -----#\n"+
            "Number of ensemble : %.0f" % (self.N_)+"\n"+
            "Number of parameter : %.0f" % (self.m_)+"\n"+
            "Number of measurements : %.0f" % (self.p_)+"\n"+
            "Standard Deviation of observations : %.6f" % np.sqrt(np.var(observations_vector,dtype=np.float64)) +"\n"+
            "#----------------------------------#\n")

        # Thinking about multiple parameters estimation
        #assert (self.p_ == len(vector_of_cell_observation)), "Dimension mismatch between your vector of grid number and your number of observations"


        # TODO sampling function to select grid parameter cell for observation


        (self.Cd_, self.D_)  = Measurements(measurements_vector,
                          self.measurements_error_percent_,
                          self.N_,
                          alpha=self.alpha_,
                          show_stats=show_parameters,
                          show_cov_matrix=show_parameters)

        """
        H = Observation_Matrix(vector_of_cell_observation,
                               number_of_grid_cell,
                               show_parameters)

        # hE Predicted value in the model observed at the same localisation as the direct measurement

         self.hE_ = H.T @ self.Ef_ # In this case the observed parameter are directly coming from my model, no forward
                                  # at this stage
        """
        self.hE_ = observations_vector



    def analysis_step(self):

        mu = np.mean(self.Ef_,1) # Computation of the ensemble mean

        A  = self.Ef_ - mu # Computation of the ensemble anomaly,
                           # individual deviation from the mean in each cell of each simulation

        hx = np.repeat(np.matrix(np.mean(self.hE_,1)).T,self.N_,axis=1) # Computation of the observed data mean


        Y  = self.hE_- hx # Computation of the observed anomaly

        import numpy.linalg as nla

        def mrdiv(b,A):
          return nla.solve(A.T,b.T).T

        C_predamp = ( Y.T @ Y)


        C  =  ( Y.T @ Y +( np.eye(self.N_)*(self.N_-1))*np.var(self.y_)*self.alpha_ ) *1.01 #inflation factor

        np.save("invert_matrix",C)

        YC = mrdiv(Y,C)

        if self.display_ :
            plt.figure(figsize=(12,8))
            plt.imshow(Y,aspect='auto')
            plt.colorbar()
            plt.title('Y')
            plt.show()

            plt.figure(figsize=(12,8))
            plt.imshow(C_predamp)
            plt.colorbar()
            plt.title('C predamp')
            plt.show()

            plt.figure(figsize=(12,8))
            plt.imshow(C)
            plt.colorbar()
            plt.title('C')
            plt.show()

            plt.figure(figsize=(12,8))
            plt.imshow(YC,aspect='auto')
            plt.colorbar()
            plt.title('YC')
            plt.show()


        KG = A @ YC.T

        self.dE = (KG @ ( self.D_ - self.hE_ ))

        self.Ea_   = self.Ef_ + self.dE

        if self.display_ :

            plt.figure(figsize=(12,8))
            plt.imshow(KG,aspect='auto')
            plt.colorbar()
            plt.title('KG')
            plt.show()


            plt.figure(figsize=(12,12))
            plt.imshow(self.dE,aspect='auto')
            plt.colorbar()
            plt.title('dE')
            plt.show()

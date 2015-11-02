'''
Created on 06.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np

class BoundedMatrix(object):
    '''
    Representation of a Matrix with finite matrix norm
    
    Attributes
    ----------
    T : numpy.ndarray([n,n], type=float)
        a numpy matrix that contains the matrix entries
    
    '''
    
    def __init__(self, T):
        self.T = T
        self.idx_computed = False
        self.evd_computed = False
        self.T_computed = True
        self.stat_idx = None
        self.normalization = ''
        
        self._compute_LDR()        
        self._compute_stationary_idx()
                
        if self.has_unique_equilibrium():
            self.set_norm('right_equilibrium')
        else:
            self.set_norm('right_2norm')
                    
        return
    
    def _compute_stationary_idx(self):
        '''
        '''
        
        if not self.idx_computed:
            # find eigenvalues to the eigenvalue of one        
                        
            self.stat_idx = np.where(np.isclose(self.eigenvalues,1))[0]
            self.idx_computed = True

    def _compute_LDR(self):
        '''
        Computes the LDR Decomposition and the stationary vectors of a matrix.
        There are two ways: Use EVD twice or only once and invert the matrix of
        '''
        
        if not self.evd_computed:
            eValues, rightEigenvectors = np.linalg.eig(self.T)            
            perm = np.argsort(eValues)[::-1]
            rightEigenvectors=rightEigenvectors[:,perm]

            eValues = eValues[perm]
            
            eValuesR, leftEigenvectors = np.linalg.eig(np.transpose(self.T))    
            permR = np.argsort(eValuesR)[::-1]
            leftEigenvectors=leftEigenvectors[:,permR].T
            
            rank = np.linalg.matrix_rank(self.T)
                        
            if rank < self.size:
                leftEigenvectors = leftEigenvectors[0:rank]
                rightEigenvectors = rightEigenvectors[:,0:rank]
                eValues = eValues[0:rank]
            
            # Check if eigenvalues agree, if not use inversion, which might be slower for sparse matrices?
            # For now we stick with good luck :)
            
            if False:
                # no more permutation needed in this case
                leftEigenvectors = np.linalg.inv(rightEigenvectors)    
                                            
            self.evals = eValues
            self.L = leftEigenvectors
            self.R = rightEigenvectors
            self.evd_computed = True
            
            self.left_norm = None
            self.right_norm = None 
        
        return
            
            
    def _compute_T(self):
        if not self.T_computed:
            if self.evd_computed:
                self.T = np.dot(np.dot(self.R, np.diagflat(self.evals)), self.L)
                self.T_computed = True
       
        return
       
    def is_normalized(self):
        return ([p for p in self._get_ev_product() if p != 1.0]) > 0
       
    @property     
    def non_zero(self):
        return 0

    @property
    def size(self):
        return self.T.shape[0]
    
    @property
    def rank(self):
        return len(self.evals)
            
    def _get_ev_product(self):
        return [np.dot(self.L[i], self.R[:,i]) for i in range( self.rank ) ]
    
    def _get_ev_left_product(self, vec = None):
        if vec is None:
            return [np.dot(self.L[i], self.L[i]) for i in range( self.rank ) ]
        else:
            return [np.dot(self.L[i]* vec, self.L[i]) for i in range( self.rank ) ]

    def _get_ev_right_product(self, vec = None):
        if vec is None:
            return [np.dot(self.R[:,i], self.R[:,i]) for i in range( self.rank ) ]
        else:
            return [np.dot(self.R[:,i] * vec, self.R[:,i]) for i in range( self.rank ) ]
                
    @property
    def stationary_idx(self):
        if not self.idx_computed:
            self._compute_stationary_idx()
            
        return self.stat_idx
        
    def set_scalar_product(self, prod):
        self.prod = prod
        
    def set_norm(self, mode):
        if mode != self.normalization:
            norm = self._get_ev_product()
            eq = self.equilibrium         
            
            
            if mode == 'left_equilibrium':
                if eq is not None:
                    norm_eq = np.sqrt(self._get_ev_left_product(1.0 / eq))
                    
                    for i in self.stationary_idx:
                        if np.sum(self.L[i]) < 0.0:
                            norm_eq[i] = -norm_eq[i]
                    self.left_norm = 1.0 / norm_eq
                    self.right_norm = 1.0 / (norm * self.left_norm)

            if mode == 'right_equilibrium':
                if eq is not None:
                    norm_eq = np.sqrt(self._get_ev_right_product(eq))

                    for i in self.stationary_idx:
                        if np.sum(self.R[:,i]) < 0.0:
                            norm_eq[i] = -norm_eq[i]
                                        
                    self.right_norm = 1.0 / norm_eq
                    self.left_norm = 1.0 / (norm * self.right_norm)

            if mode == 'left_2norm':
                norm_eq = np.sqrt(self._get_ev_left_product())
                for i in self.stationary_idx:
                    if np.sum(self.L[i]) < 0.0:
                        norm_eq[i] = -norm_eq[i]
                
                    self.left_norm = 1.0 / norm_eq
                    self.right_norm = 1.0 / (norm * self.left_norm)

            if mode == 'right_2norm':
                norm_eq = np.sqrt(self._get_ev_right_product())
                
                for i in self.stationary_idx:
                    if np.sum(self.R[:,i]) < 0.0:
                        norm_eq[i] = -norm_eq[i]
                
                self.right_norm = 1.0 / norm_eq
                self.left_norm = 1.0 / (norm * self.right_norm)

            self.normalization = mode
                          
        return
        
    @property
    def eigenvalues(self):
        return self.evals

    @property
    def left_eigenvectors(self):
        if self.left_norm is not None:
            return np.dot( np.diagflat(self.left_norm), self.L)
        else:
            return self.L

    @property
    def right_eigenvectors(self):
        if self.right_norm is not None:
            return np.dot(self.R, np.diagflat(self.right_norm))
        else:
            return self.R
    
    def has_unique_equilibrium(self):
        '''
        '''
        return len(self.stationary_idx) == 1

    @property
    def equilibrium(self):
        '''
        '''
        
        if self.has_unique_equilibrium():
            eq = self.L[self.stationary_idx[0]]
            
            eq = eq / np.sum(eq)
            return eq
        
        return None
    
    def range_sum(self, r):
        beg = r[0]
        end = r[1]
        
        if type(end) is int:
            arr = [ (np.power(l, end + 1) - np.power(l, beg)) / (l - 1.0) if l < 1.0 else end - beg + 1 for l in self.evals ]
            di = np.array( arr )
        else:
            di = np.zeros(self.size)
        return np.dot(np.dot(self.right_eigenvectors, np.diagflat(di)), self.left_eigenvectors)

    
    @property
    def left_steady_vectors(self):
        return self.left_eigenvectors[self.stationary_idx]

    def right_steady_vectors(self):
        return self.right_eigenvectors[self.stationary_idx]
    
    @property
    def left_steady(self):
        return self.lvecs
    
                
class MSM(BoundedMatrix):
    '''
    Representation of a Markov State Model. Allows for several ways to compute properties of a MSM
    
    Parameters
    T : BoundedMatrix
    '''


    def __init__(self, T):
        super(MSM, self).__init__(T)
        
    def __repr__(self, *args, **kwargs):
        return object.__repr__(self, *args, **kwargs)
    
    def __str__(self, *args, **kwargs):
        return object.__str__(self, *args, **kwargs)
        
    def _tau(self, idx):
        if idx =='all':
            return self.T
        else:
            return np.dot(np.diagflat(self._sv(idx)), self.T)
    
    def _sv(self, c):
        s = np.zeros(self.size)
        
        if c is not list:
            c = [c]
        
        for i in c:
            s[i] = 1
        return s

    def path_probabilities(self, description):
#        s = self._sv(description[0])
        s = np.ones(self.size)
        
        for i in description:
            if type(i) is tuple:
                if type(i[1]) is int:
                    s = np.dot(s, np.linalg.matrix_power(self._tau(i[0]), i[1]))
                if type(i[1]) is str:
                    s = np.dot(s, np.linalg.matrix_power(self._tau(i[0]), 1))
                elif type(i[1]) is tuple:
                    s = np.dot(s, BoundedMatrix(self._tau(i[0])).range_sum(i[1]))
            else:
                s = np.dot(s, self._tau(i))

        return np.dot(s, np.ones(self.size))
    
    def OOM_representation(self):
        return
    
    def meanfirstpassagetimes(self, target):
        return
    
    def mpft_matrix(self):
        return
    
    def committor(self):
        return
    
    def nearest_reversible(self):
        return
    
    def new_state_discover(self):
        return
    
    def reactive_flux(self):
        return
    
    def net_flux(self):
        return
    
    def reactivity(self):
        return
    
    def state_exit_flow_rate_sensitivity(self):
        return
    
    def equilibrium_sensitivity(self):
        return
    
    def steady_state_equilibrium(self):
        return
    
    def implied_timescales(self):
        return
    
    def NC_LDR_decomposition(self):
        return
    
    def transition_rate(self):
        return
    
    def probability_current(self):
        return
    
    def forward_committor(self):
        return
    
    def backward_committor(self):
        return
    
    def _matrix_range_sum(self, range):
        return
    
    def _matrix_range_mean_sum(self, range):
        return
    
    def _steady_eigenvector(self):
        return
    
    def _to_LDR(self):
        return
    
    def _from_LDR(self):
        return
    
    def correlation(self, obd1, obs2):
        return
    
    def autocorrelation(self, obs):
        return
    
    def is_stochastic(self):
        return
    
    def is_symmetric(self):
        return
    
    def has_row_sum_one(self):
        return
    
    def has_row_sum_zero(self):
        return
    
    def has_total_sum_one(self):
        return
    
    def spectral_matrix_norm(self):
        return
    
    def frobenius_norm(self):
        return
    
    def matrix_norm(self):
        return
    
    def normalize_LDR(self):
        return
    
    def as_OOM(self):
        return
    
    def as_PMM(self):
        return
    
    def get_TM(self):
        return
    
    def get_XM(self):
        return
    
    def get_CM(self):
        return
            
    def limited_committor(self):
        return
    
    def limited_mfpt(self):
        return
    
    def forward_committor_sensitivity(self, A, B, index):
        """ 
        calculate the sensitivity matrix for index of the forward committor from A to B given transition matrix T.
        
        Parameters
        ----------
        T : np.ndarray shape = (n, n)
            Transition matrix
        A : array like
            List of integer state labels for set A
        B : array like
            List of integer state labels for set B
        index : int     
            entry of the committor for which the sensitivity is to be computed

        Returns
        -------
        x : ndarray, shape=(n, n)
            Sensitivity matrix for entry index around transition matrix T. Reversibility is not assumed.
        """
        
        n = len(self.T)
        set_X = np.arange(n)#set(range(n))
        set_A = np.unique(A)#set(A)
        set_B = np.unique(B)#set(B)
        set_AB = np.union1d(set_A, set_B)#set_A | set_B
        notAB = np.setdiff1d(set_X, set_AB, True)#list(set_X - set_AB)
        m = len(notAB)
    
        K = self.T - np.diag(np.ones(n))
    
        U = K[np.ix_(notAB, notAB)]
        
        v = np.zeros(m)
        
        #for i in xrange(0, m):
        #   for k in xrange(0, len(set_B)):
        #       v[i] = v[i] - K[notAB[i], B[k]]
        v[:] = v[:] - K[notAB[:], B[:]]
    
        qI = np.linalg.solve(U, v)
    
        q_forward = np.zeros(n)
        #q_forward[set_A] = 0 # double assignment.
        q_forward[set_B] = 1
        #for i in range(len(notAB)):
        q_forward[notAB[:]] = qI[:]
            
        target = np.eye(1,n,index)
        target = target[0,notAB]
    
        UinvVec = np.linalg.solve(U.T, target)
        Siab = np.zeros((n,n))
            
        for i in xrange(m):
            Siab[notAB[i]] = - UinvVec[i] * q_forward
    
        return Siab
    
    def backward_committor_sensitivity(self, A, B, index):
        """ 
        calculate the sensitivity matrix for index of the backward committor from A to B given transition matrix T.
        
        Parameters
        ----------
        T : np.ndarray shape = (n, n)
            Transition matrix
        A : array like
            List of integer state labels for set A
        B : array like
            List of integer state labels for set B
        index : int 
            entry of the committor for which the sensitivity is to be computed
            
        Returns
        -------
        x : ndarray, shape=(n, n)
            Sensitivity matrix for entry index around transition matrix T. Reversibility is not assumed.
            
        """
        
        # This is really ugly to compute. The problem is, that changes in T induce changes in
        # the stationary distribution and so we need to add this influence, too
        # I implemented something which is correct, but don't ask me about the derivation
        
        n = len(self.T)
        
        trT = np.transpose(self.T)
        
        one = np.ones(n)
        eq = self.stationary_distribution(self.T)
        
        mEQ = np.diag(eq)
        mIEQ = np.diag(1.0 / eq)
        mSEQ = np.diag(1.0 / eq / eq)
        
        backT = np.dot(mIEQ, np.dot( trT, mEQ))
        
        qMat = self.forward_committor_sensitivity(A, B, index)
        
        matA = trT - np.identity(n)
        matA = np.concatenate((matA, [one]))
        
        phiM = np.linalg.pinv(matA)
        
        phiM = phiM[:,0:n]
        
        trQMat = np.transpose(qMat)
        
        d1 = np.dot( mSEQ, np.diagonal(np.dot( np.dot(trT, mEQ), trQMat), 0) )
        d2 = np.diagonal(np.dot( np.dot(trQMat, mIEQ), trT), 0)
            
        psi1 = np.dot(d1, phiM)
        psi2 = np.dot(-d2, phiM)
        
        v1 = psi1 - one * np.dot(psi1, eq)
        v3 = psi2 - one * np.dot(psi2, eq)
        
        part1 = np.outer(eq, v1)
        part2 = np.dot( np.dot(mEQ, trQMat), mIEQ)
        part3 = np.outer(eq, v3)
        
        sensitivity = part1 + part2 + part3
        
        return sensitivity
    
    def eigenvalue_sensitivity(self, k):
        """ 
        calculate the sensitivity matrix for eigenvalue k given transition matrix T.
        
        Parameters
        ----------
        T : np.ndarray shape = (n, n)
            Transition matrix
        k : int
            eigenvalue index for eigenvalues order descending
            
        Returns
        -------
        x : ndarray, shape=(n, n)
            Sensitivity matrix for entry index around transition matrix T. Reversibility is not assumed.
            
        """
            
        eValues, rightEigenvectors = np.linalg.eig(self.T)
        leftEigenvectors = np.linalg.inv(rightEigenvectors)    
        
        perm = np.argsort(eValues)[::-1]
    
        rightEigenvectors=rightEigenvectors[:,perm]
        leftEigenvectors=leftEigenvectors[perm]
        
        sensitivity = np.outer(leftEigenvectors[k], rightEigenvectors[:,k])
        
        return sensitivity
    
    def timescale_sensitivity(self, k):
        """ 
        calculate the sensitivity matrix for timescale k given transition matrix T.
        
        Parameters
        ----------
        T : np.ndarray shape = (n, n)
            Transition matrix
        k : int
            timescale index for timescales of descending order (k = 0 for the infinite one)
            
        Returns
        -------
        x : ndarray, shape=(n, n)
            Sensitivity matrix for entry index around transition matrix T. Reversibility is not assumed.
        """
            
        eValues, rightEigenvectors = np.linalg.eig(self.T)
        leftEigenvectors = np.linalg.inv(rightEigenvectors)    
        
        perm = np.argsort(eValues)[::-1]
    
        eValues = eValues[perm]
        rightEigenvectors=rightEigenvectors[:,perm]
        leftEigenvectors=leftEigenvectors[perm]
        
        eVal = eValues[k]
        
        sensitivity = np.outer(leftEigenvectors[k], rightEigenvectors[:,k])
        
        if eVal < 1.0:
            factor = 1.0 / (np.log(eVal)**2) / eVal
        else:
            factor = 0.0    
        
        sensitivity *= factor
        
        return sensitivity
    
    # TODO: The eigenvector sensitivity depends on the normalization, e.g. l^T r = 1 or norm(r) = 1
    # Should we fix that or add another option. Also the sensitivity depends on the initial eigenvectors
    # Now everything is set to use norm(v) = 1 for left and right
    # In the case of the stationary distribution we want sum(pi) = 1, so this function
    # does NOT return the same as stationary_distribution_sensitivity if we choose k = 0 and right = False!
    
    # TODO: If we choose k = 0 and right = False we might throw a warning!?!
    
    def eigenvector_sensitivity(self, k, j, right=True):
        """ 
        calculate the sensitivity matrix for entry j of left or right eigenvector k given transition matrix T.
        
        Parameters
        ----------
        T : np.ndarray shape = (n, n)
            Transition matrix
        k : int
            eigenvector index ordered with descending eigenvalues
        j : int
            entry of eigenvector k for which the sensitivity is to be computed
        right : boolean (default: True)
            If set to True (default) the right eigenvectors are considered, otherwise the left ones
            
        Returns
        -------
        x : ndarray, shape=(n, n)
            Sensitivity matrix for entry index around transition matrix T. Reversibility is not assumed.
        
        Notes
        -----
        Eigenvectors can naturally be scaled and so will their sensitivity depend on their size.
        For that reason we need to agree on a normalization for which the sensitivity is computed.
        Here we use the natural norm(vector) = 1 condition which is different from the results, e.g.
        the rdl_decomposition returns. 
        This is especially important for the stationary distribution, which is the first left eigenvector.
        For this reason this function return a different sensitivity for the first left eigenvector
        than the function stationary_distribution_sensitivity and this function should not be used in this
        case!
        
        """
        
        n = len(self.T)
        
        if not right:    
            T = np.transpose(self.T)
        
        eValues, rightEigenvectors = np.linalg.eig(self.T)
        leftEigenvectors = np.linalg.inv(rightEigenvectors)        
        perm = np.argsort(eValues)[::-1]
    
        eValues = eValues[perm]
        rightEigenvectors=rightEigenvectors[:,perm]
        leftEigenvectors=leftEigenvectors[perm]
            
        rEV = rightEigenvectors[:,k]
        lEV = leftEigenvectors[k]
        eVal = eValues[k]
        
        vecA = np.zeros(n)
        vecA[j] = 1.0
               
        matA = T - eVal * np.identity(n)
            # Use here rEV as additional condition, means that we assume the vector to be
            # orthogonal to rEV
        matA = np.concatenate((matA, [rEV]))
                    
        phi = np.linalg.lstsq(np.transpose(matA), vecA)    
            
        phi = np.delete(phi[0], -1)
                    
        sensitivity = -np.outer(phi,rEV) + np.dot(phi,rEV) * np.outer(lEV, rEV) 
        
        if not right:
            sensitivity = np.transpose(sensitivity)          
            
        return sensitivity
    
    def stationary_distribution_sensitivity(self, j):
        r"""Calculate the sensitivity matrix for entry j the stationary distribution vector given transition matrix T.
    
        Parameters
        ----------
        T : np.ndarray shape = (n, n)
            Transition matrix
        j : int
            entry of stationary distribution for which the sensitivity is to be computed
            
        Returns
        -------
        x : ndarray, shape=(n, n)
            Sensitivity matrix for entry index around transition matrix T. Reversibility is not assumed.
        
        Notes
        -----
        Note, that this function uses a different normalization convention for the sensitivity compared to
        eigenvector_sensitivity. See there for further information.
        """
                
        n = len(self.T)
            
        lEV = np.ones(n)
        rEV = self.stationary_distribution(self.T)
        eVal = 1.0
        
        T = np.transpose(self.T)
        
        vecA = np.zeros(n)
        vecA[j] = 1.0
                   
        matA = T - eVal * np.identity(n)
        # normalize s.t. sum is one using rEV which is constant
        matA = np.concatenate((matA, [lEV]))
                    
        phi = np.linalg.lstsq(np.transpose(matA), vecA)    
        phi = np.delete(phi[0], -1)
                        
        sensitivity = -np.outer(rEV, phi) + np.dot(phi,rEV) * np.outer(rEV, lEV)           
            
        return sensitivity
    
    def mfpt_sensitivity(self, target, j):
        """ 
        calculate the sensitivity matrix for entry j of the mean first passage time (MFPT) given transition matrix T.
        
        Parameters
        ----------
        T : np.ndarray shape = (n, n)
            Transition matrix
        target : int
            target state to which the MFPT is computed
        j : int
            entry of the mfpt vector for which the sensitivity is to be computed
            
        Returns
        -------
        x : ndarray, shape=(n, n)
            Sensitivity matrix for entry index around transition matrix T. Reversibility is not assumed.
        """
        
        n = len(self.T)
        
        matA = self.T - np.diag(np.ones((n)))
        matA[target] *= 0
        matA[target, target] = 1.0
        
        tVec = -1. * np.ones(n);
        tVec[target] = 0;
        
        mfpt = np.linalg.solve(matA, tVec)
        aVec = np.zeros(n)
        aVec[j] = 1.0
        
        phiVec = np.linalg.solve(np.transpose(matA), aVec )
        
        # TODO: Check sign of sensitivity!
            
        sensitivity = -1.0 * np.outer(phiVec, mfpt)
        sensitivity[target] *= 0;
        
        return sensitivity
    
    def expectation_sensitivity(self, a):
        r"""Sensitivity of expectation value of observable A=(a_i).
    
        Parameters
        ----------
        T : (M, M) ndarray
            Transition matrix
        a : (M,) ndarray
            Observable, a[i] is the value of the observable at state i.
    
        Returns
        -------
        S : (M, M) ndarray
            Sensitivity matrix of the expectation value.
        
        """
        n=self.T.shape[0]
        S=np.zeros((n, n))
        for i in range(n):
            S+=a[i]*self.stationary_distribution_sensitivity(i)
        return S
    
    @property
    def all(self):
        return range(self.size)
        
X = MSM(np.array([[0.9,0.02,0.08], [0.03,0.92,0.05], [0.01,0.02,0.97]]))

#X._compute_LDR()

#print X.left_eigenvectors    
#print X.right_eigenvectors
#X.set_norm('left_2norm')
#print X.left_eigenvectors
#print X.right_eigenvectors

#print np.dot(X.left_eigenvectors, X.right_eigenvectors)
#X._compute_T()

#print X.T

#print X.forward_committor_sensitivity([0], [2], 1)

All = 'all'

print X.path_probabilities([0,([0,1], (0,1)),0])
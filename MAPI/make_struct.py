import numpy as np

class distorted_crystal:
    
    def __init__(self):
        
        self.oct_1 = np.array([[0.5,0.5,0.0]
                               [0.0,0.0,0.0],
                               [1.0,0.0,0.0],
                               [0.0,1.0,0.0],
                               [1.0,1.0,0.0]])
        self.oct_2 = np.array([[0.5,0.5,0.0]]) # only 1 atom that isn't connected by PBC

    def make_supercell(self,num_layers=3,phis=[25,-3,-15]):
        
        pass

    def _rotate_layer(self,phi):

        rad = phi*np.pi/180 # go to radians

        # ccw rotation matrix (in x-y plane only)
        R = np.array([[ np.cos(rad),-np.sin(rad), 0.0],
                      [ np.sin(rad), np.cos(rad), 0.0],
                      [         0.0,         0.0, 1.0]])

        _oct_1 = np.zeros(self.oct_1.shape)
        _oct_2 = np.zeros(self.oct_2.shape)

        for ii in range(_oct_1.shape[0]):
            _oct_1[ii,:] = np.matmul(R,self.oct_1[ii,:])
        for ii in range(_oct_1.shape[0]):
            _oct_1[ii,:] = np.matmul(R,self.oct_1[ii,:])


        




import numpy as np 

import matplotlib.pyplot as plt


class c_bragg_scattering:

    def __init__(self,basis_pos,scattering_lens,num_angles,reps=None):

        """
        calculates bragg scattering intensity as a function of scattering angle in the xy
        plane
        """

        self.basis_pos = np.atleast_2d(np.array(basis_pos))
        self.num_basis = self.basis_pos.shape[0]

        self.num_angles = num_angles
        self.angles = np.linspace(0,2*np.pi-(2/np.pi/num_angles),num_angles)

        self.scattering_lens = np.array(scattering_lens)

        if reps is not None:
            
            self.finite = True
            self.reps = np.array(reps)
            self.num_reps = np.prod(self.reps)

            _x, _y = np.meshgrid(np.arange(0,reps[0]),np.arange(0,reps[1]),indexing='ij')
            self.unitcell_coords = np.c_[_x.flatten(),_y.flatten()].astype(float)
 
        else:
            self.finite = False
    
    def calculate_scattering(self,d_spacing=1):

        """
        Q = k-k' = k [x - x'] = k[(1-cos(theta))x + sin(theta)y]

        scattering only occurs for k = 2 pi / d
        
        with d the "d spacing", i.e. the distance between a set of lattice planes. 
        d_1 = R, 
        """

        self.k = 2*np.pi/d_spacing

        self.Q = np.zeros((self.num_angles,2),dtype=float)
        self.Q[:,0] = self.k*(1-np.cos(self.angles))
        self.Q[:,1] = -self.k*np.sin(self.angles)

        self.structure_factor = np.zeros(self.num_angles)
        self.form_factor = np.zeros(self.num_angles)
        self.intensity = np.zeros(self.num_angles)

        for ii in range(self.num_angles):
            
            print(ii)

            _Q = np.tile(self.Q[ii].reshape(1,2),reps=(self.num_basis,1))
            _exp_Qt = np.exp(1j*np.sum(_Q*self.basis_pos,axis=1))
            _ffac = np.sum(self.scattering_lens*_exp_Qt)

            _Q = np.tile(self.Q[ii].reshape(1,2),reps=(self.num_reps,1))
            _exp_QR = np.exp(1j*np.sum(_Q*self.unitcell_coords,axis=1))
            _sfac = np.sum(_exp_QR)

            self.structure_factor[ii] = np.abs(_sfac)**2
            self.form_factor[ii] = np.abs(_ffac)**2
            self.intensity[ii] = np.abs(_sfac*_ffac)**2



scatt = c_bragg_scattering(basis_pos=[0,0],scattering_lens=[1.0],num_angles=3599,reps=[100,100])
scatt.calculate_scattering(d_spacing=2)#/np.sqrt(2))#1.09)


angles = scatt.angles #* 180 / np.pi
intensity = scatt.intensity
form_factor = scatt.form_factor
structure_factor = scatt.structure_factor

fig, ax = plt.subplots()#subplot_kw={'projection':'polar'})

plt.plot(angles/np.pi,intensity,c='k')
plt.plot(angles/np.pi,form_factor,c='b')
plt.plot(angles/np.pi,structure_factor,c='m')
#plt.plot(2
plt.show()


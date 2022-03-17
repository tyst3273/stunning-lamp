import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class distorted_crystal:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self):

        self.a = 6.002
        self.basis = np.array([[0.25,0.75,0.50],  # Cs1
                               [0.25,0.25,0.00],  # Pb1 
                               [0.00,0.00,0.00],  # Br11
                               [0.50,0.00,0.00],  # Br12
                               [0.00,0.50,0.00],  # Br13
                               [0.25,0.25,0.50],  # Br14
                               [0.50,0.50,0.00],  # Br15
                               [0.75,0.25,0.50],  # Cs2
                               [0.75,0.75,0.00],  # Pb2
                               [0.75,0.75,0.50]]) # Br21
        self.basis_types = ['Cs', 'Pb','Br','Br', 'Br','Br','Br', 'Cs', 'Pb','Br']
        self.oct_1 = [1,2,3,4,5,6]
        self.oct_2 = [8,9]
        self.n_atom = self.basis.shape[0]

    # ----------------------------------------------------------------------------------------------

    def make_supercell(self,ang_deg=[25,-3,-15]):

        self.ang_deg = ang_deg
        num_layers = len(ang_deg)
        self.num_layers = num_layers

        _pos = np.zeros((self.n_atom*num_layers,3)) 
        self.types = np.zeros(self.n_atom*num_layers).astype(str)

        # rotate and stack layers
        work = np.copy(self.basis)
        for ii in range(num_layers):
            shift = ii*self.n_atom

            # rotate first octahedra
            _tmp = work[self.oct_1,:]
            _tmp = self._rotate(_tmp,ang_deg[ii])
            work[self.oct_1,:] = _tmp[:,:]

            # rotate 2nd one. first, shift it to origin.
            _tmp = work[self.oct_2,:]
            _tmp[:,:2] = _tmp[:,:2]-0.5
            _tmp = self._rotate(_tmp,-ang_deg[ii]) # now rotate
            # then shift it back so that it shares corner with 1st octahedra
            _tmp[:,0] = _tmp[:,0]+work[self.oct_1[-1],0]
            _tmp[:,1] = _tmp[:,1]+work[self.oct_1[-1],1]
            work[self.oct_2,:] = _tmp[:,:]

            # now shift along z axis
            work[:,2] = work[:,2]+ii
            work[:,0] = work[:,0]-work[:,0].min()
            work[:,1] = work[:,1]-work[:,1].min()

            _pos[shift:shift+self.n_atom,:] = work[:,:]
            self.types[shift:shift+self.n_atom] = self.basis_types[:]

            # now reset positions and do it again
            work[:,:] = self.basis[:,:]

        # now normalize so that it all is in crystal coords
        _pos[:,2] = _pos[:,2]/self.num_layers

        sort = np.argsort(self.types)
        cs = np.argwhere(self.types == 'Cs').flatten()
        self.n_cs = cs.size
        pb = np.argwhere(self.types == 'Pb').flatten()
        self.n_pb = pb.size
        br = np.argwhere(self.types == 'Br').flatten()
        self.n_br = br.size
        self.pos = np.zeros((self.n_atom*num_layers,3))
        self.pos[:self.n_cs,:] = _pos[cs,:]
        self.pos[self.n_cs:self.n_cs+self.n_pb,:] = _pos[pb,:]
        self.pos[self.n_cs+self.n_pb:self.n_cs+self.n_pb+self.n_br,:] = _pos[br,:]
            
    # ----------------------------------------------------------------------------------------------

    def write_poscar(self,outfile='POSCAR'):
        with open(outfile,'w') as fout:
            _ = 'ANGLES(degrees) '
            for ii in range(self.num_layers):
                _ = _ + f' {self.ang_deg[ii]:4.2f}'
            fout.write(_+'\n')
            fout.write(f' {self.a:8.6f}\n')
            fout.write(f' {np.sqrt(2):6.5f}  0.00000  0.00000\n')
            fout.write(f' 0.00000  {np.sqrt(2):6.5f}  0.00000\n')
            fout.write(f' 0.00000  0.00000  {self.num_layers:6.5f}\n')
            fout.write('Cs Pb Br\n')
            fout.write(f' {self.n_cs} {self.n_pb} {self.n_br}\n')
            fout.write('Direct\n')
            for ii in range(self.pos.shape[0]):
                fout.write(f' {self.pos[ii,0]:9.8f} {self.pos[ii,1]:9.8f} {self.pos[ii,2]:9.8f}\n')

    # ----------------------------------------------------------------------------------------------

    def _rotate(self,_in,phi):

        rad = phi*np.pi/180 # go to radians

        # ccw rotation matrix (in x-y plane only)
        R = np.array([[ np.cos(rad),-np.sin(rad), 0.0],
                      [ np.sin(rad), np.cos(rad), 0.0],
                      [         0.0,         0.0, 1.0]])

        _out = np.zeros(_in.shape)
        for ii in range(_in.shape[0]):
            _out[ii,:] = np.matmul(R,_in[ii,:])

        return _out

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    crystal = distorted_crystal()

    crystal.make_supercell(ang_deg=[-12,0,12])
    crystal.write_poscar('POSCAR_1')

    crystal.make_supercell(ang_deg=[-12,12,-10])
    crystal.write_poscar('POSCAR_2')

    crystal.make_supercell(ang_deg=[-8,-8,10])
    crystal.write_poscar('POSCAR_3')

    crystal.make_supercell(ang_deg=[4,12,16])
    crystal.write_poscar('POSCAR_4')






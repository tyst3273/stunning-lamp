import numpy as np
from scipy.linalg import eigh

def solve(nb=2):

    nk = 100
    k = np.linspace(-1/2,1/2,nk,endpoint=False)

    nx = 101
    x = np.linspace(0,1,nx,endpoint=False)

    V0 = 5
    V = V0*np.cos(2*np.pi*x)/2

    VG = np.fft.fftfreq(nx,x[1]-x[0])
    V_fft = np.fft.fft(V,norm='forward').round(6)

    Gcut = VG.max()/2
    G = VG[np.flatnonzero(np.abs(VG) <= Gcut)]
    nG = G.size

    KE = np.zeros((nG,nG))
    PE = np.zeros((nG,nG),dtype=complex)

    ek = np.zeros((nk,nb))
    wk = np.zeros((nk,nb,nG),dtype=complex)

    ek0 = np.zeros((nk,nb))
    wk0 = np.zeros((nk,nb,nG),dtype=complex)

    for kk in range(nk):

        KE[...] = 0.0
        PE[...] = 0.0

        for ii in range(nG):
            for jj in range(ii,nG):

                if ii == jj:
                    KE[ii,jj] = 4*np.pi**2*0.5*(k[kk]+G[ii])**2

                Gii = G[ii]; Gjj = G[jj]
                GG = Gii-Gjj
                ind = np.flatnonzero(VG == GG)[0]
                PE[ii,jj] = V_fft[ind]

        evals, evecs = eigh(KE,subset_by_index=[0,nb-1])
        ek0[kk,:] = evals
        wk0[kk,...] = evecs.T

        evals, evecs = eigh(KE+PE,subset_by_index=[0,nb-1],lower=False)
        ek[kk,:] = evals
        wk[kk,...] = evecs.T

    return ek0, wk0, ek, wk, x, k, G, V, V0


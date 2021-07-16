from scipy import *
import scipy.weave as weave
from pylab import *
import sys #, symeig
import excor
from scipy import linalg
from scipy import optimize
from scipy import integrate
from scipy import interpolate
from scipy import special
from scipy import misc

##############################################################
###     Everything connected with fcc crystal structure  #####
###     Could be generalized to an arbitrary crystal    #####
##############################################################
class FccLattice:
    "Handles fcc crystal. Could be generalized to arbitrary crystal."
    def __init__(self, LatConst):
        self.a0 = array([0.5*LatConst,0.5*LatConst,0])
        self.a1 = array([0.5*LatConst,0,0.5*LatConst])
        self.a2 = array([0,0.5*LatConst,0.5*LatConst])
        Vc = dot(cross(self.a0,self.a1),self.a2) # Volume
        self.Volume = abs(Vc)
        print "Volume is ", self.Volume
        self.b0 = (2*pi/Vc)*cross(self.a1,self.a2)
        self.b1 = (2*pi/Vc)*cross(self.a2,self.a0)
        self.b2 = (2*pi/Vc)*cross(self.a0,self.a1)
        # Special points in Brillouin zone
        brs = 2*pi/LatConst
        self.GPoint = [0,0,0]
        self.LPoint = array([0.5,  0.5,  0.5])*brs
        self.KPoint = array([0.75, 0.75, 0])*brs
        self.XPoint = array([1.0,  0.0,  0])*brs
        self.WPoint = array([1.0,  0.5,  0])*brs
        
    def RMuffinTin(self):
        return 0.5*sqrt(dot(self.a0,self.a0)) # Spheres just touch

    def GenerateReciprocalVectors(self, q, CutOffK):
        # Many reciprocal vectors are generated and later only the shortest are used
        Kmesh0=[]
        for n in range(-q,q+1):
            for l in range(-q,q+1):
                for m in range(-q, q+1):
                    vec = n*self.b0+l*self.b1+m*self.b2
                    if dot(vec,vec) <= CutOffK**2:
                        Kmesh0.append(vec)
                    
        Kmesh0.sort(lambda x,y: cmp(dot(x,x), dot(y,y)))
        self.Km = array(Kmesh0)
        print "K-mesh size=", len(self.Km)
        
    def ChoosePointsInFBZ(self, nkp, type=0): # Chooses the path in the 1BZ we will use
        
        def kv0(iq, q):
            return (iq-int((q+1.5)/2)+1)/(q+0.0)
        
        if type==0: # Choose mesh in the 1BZ to cover the whole space - for SC calculation
            kp=[]
            for i0 in range(nkp):
                r0 = kv0(i0,nkp)
                print 'r0=', r0
                for i1 in range(nkp):
                    r1 = kv0(i1,nkp)
                    for i2 in range(nkp):
                        r2 = kv0(i2,nkp)
                        k = self.b0*r0+self.b1*r1+self.b2*r2
                        kp.append(k)
            print "Number of all k-points=", len(kp)

            kpc = []
            for k in kp:
                kpc.append(sort(k))
                
            # ChooseIrreducible k-points only
            # The function performs all symmetry operations of a cubic point-group to each k-point and
            # keeps only thos k-points which can not be obtained from another k-point by group operation.
            # These k-points are obviously irreducible.
            irkp = []       # temporary list where irreducible k points will be stored
            wkp  = []       # temporary weights
            while len(kpc)>0: # continues until all k-points are grouped into irreducible classes 
                tk = kpc[0]               # we concentrate on the k-point which is the first in the list
                irkp.append(tk)          # the first can be stored as irreducible
                wkp.append(0)            # and the weights for this irreducible k-point is set to zero
                # We go over 48 symmetry operations of cubic system:
                # Each wector component can change sign: 2^3=8 possibilities
                # All permutations of components: 3!=6
                # Since the operations are independent, we have 3!*2^3=48 operations == number of cubic point group operations
                for ix in [-1,1]:  # three loops for all possible sign changes 
                    for iy in [-1,1]:
                        for iz in [-1,1]:
                            nk = sort([ix*tk[0], iy*tk[1], iz*tk[2]]) # sorted so that we do not need to try all permutations
                            ii=0
                            while ii<len(kpc): # This permutation and sign change leads to some element still in the list of k-points?
                                diff = sum(abs(nk - kpc[ii]))
                                if diff<1e-6:
                                    del kpc[ii] # These two k-points are the same
                                    wkp[-1] += 1.
                                else:
                                    ii+=1

            # irreducible k-points are stored in the output vectors
            self.wkp = array(wkp)/sum(wkp)
            self.kp = array(irkp)

            print "Number of irreducible k points is: ", len(self.kp)
            #for ik,k in enumerate(self.kmesh):
            #    print "%10.6f"*3 % tuple(k), '  ', self.wkp[ik]
            
        else:        # Choose one particular path in the 1BZ - for plotting purposes
            nkp = 4*int(nkp/4.)+1
            print "nkp=", nkp
            self.kp = zeros((nkp,3), dtype=float)
            N0=nkp/4

            self.Points = [('$\Gamma$', 0), ('$X$', N0), ('$L$', 2*N0), ('$\Gamma$', 3*N0), ('$K$', 4*N0)]
            for i in range(N0): self.kp[i,:]      = self.GPoint + (self.XPoint-self.GPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0+i,:]   = self.XPoint + (self.LPoint-self.XPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0*2+i,:] = self.LPoint + (self.GPoint-self.LPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0*3+i,:] = self.GPoint + (self.KPoint-self.GPoint)*i/(N0-0.)
            self.kp[4*N0] = self.KPoint
            

########################################################
# Routines for solving the ODE problem and Poisson EQ. #
########################################################
def Numerov(F, dx, f0=0.0, f1=1e-3):
    codeNumerov="""
      double h2 = dx*dx;
      double h12 = h2/12;
      
      double w0 = (1-h12*F(0))*Solution(0);
      double Fx = F(1);
      double w1 = (1-h12*Fx)*Solution(1);
      double Phi = Solution(1);
      
      double w2;
      for (int i=2; i<Nmax; i++){
        w2 = 2*w1 - w0 + h2*Phi*Fx;
        w0 = w1;
        w1 = w2;
        Fx = F(i);
        Phi = w2/(1-h12*Fx);
        Solution(i) = Phi;
      }
    """
    Nmax = len(F)
    dx = float(dx)
    Solution = zeros(Nmax, dtype=float)
    Solution[0] = f0
    Solution[1] = f1
    weave.inline(codeNumerov, ['F', 'Nmax', 'dx', 'Solution'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return Solution

def NumerovGen(F, U, dx, f0=0.0, f1=1e-3):
    codeNumerov="""
      double h2 = dx*dx;
      double h12 = h2/12;
      
      double w0 = Solution(0)*(1-h12*F(0))-h12*U(0);
      double Fx = F(1);
      double Ux = U(1);
      double w1 = Solution(1)*(1-h12*Fx)-h12*Ux;
      double Phi = Solution(1);
      
      double w2;
      for (int i=2; i<Nmax; i++){
        w2 = 2*w1 - w0 + h2*(Phi*Fx+Ux);
        w0 = w1;
        w1 = w2;
        Fx = F(i);
        Ux = U(i);
        Phi = (w2+h12*Ux)/(1-h12*Fx);
        Solution(i) = Phi;
      }
    """
    Nmax = len(F)
    dx = float(dx)
    Solution = zeros(Nmax, dtype=float)
    Solution[0] = f0
    Solution[1] = f1
    weave.inline(codeNumerov, ['F', 'U', 'Nmax', 'dx', 'Solution'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return Solution

def CRHS(E, l, R, Veff):
    "RHS for solving the Schroedinger equations by Numerov. To achive sufficient speed, uses C++."
    codeRHS="""
        for (int i=0; i<N; i++){
           RHS(i) = 2*( -E + 0.5*l*(l+1)/(R(i)*R(i)) + Veff(i) );
        }
    """
    N = len(R)
    RHS = zeros(len(R), dtype=float)
    weave.inline(codeRHS, ['N', 'E', 'l', 'R', 'Veff', 'RHS'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return RHS

def SolvePoisson(Zq, R, rho):
    """Given the input density rho, calculates the Hartree potential
    The boundary conditions used are U(0)=0 and U(S)=Zq. The boundary condition at S is only a constant shift
    of potential which is later readjusted by choosing MT zero. So, this boundary condition is choosen for convenience.
    """
    codeNumerovInh="""
        double h2 = dx*dx;
        double h12 = h2/12;
        
        double w0 = Solution(0)-h12*U(0);
        double Ux = U(1);
        double w1 = Solution(1)-h12*Ux;
        double Phi = Solution(1);
        
        double w2;
        for (int i=2; i<Nmax; i++){
          w2 = 2*w1 - w0 + h2*Ux;
          w0 = w1;
          w1 = w2;
          Ux = U(i);
          Phi = w2+h12*Ux;
          Solution(i) = Phi;
        }
     """
    U = array([-4*pi*r*rho[i] for i,r in enumerate(R)])
    Nmax = len(R)
    dx = float( (R[-1]-R[0])/(len(R)-1.) )
    Solution = zeros(len(R), dtype=float)
    
    Solution[0]=0
    Solution[1]=(R[1]-R[0]) # Boundary condition for U_H=V_H/r
    
    weave.inline(codeNumerovInh, ['U', 'Nmax', 'dx', 'Solution'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    # adding homogeneous solution to satisfay boundary conditions: U(0)=0, U(infinity)=Z
    alpha = (Zq - Solution[-1])/R[-1]
    Solution += alpha*R
    return Solution

#############################
#       LAPW Routins        #
#############################
def ComputeInterstitialOverlap(Km, RMuffinTin, Vol):
    """ Overlap in the interstitials can be calculated outside the k-loop
    
        Please see Eq.46 on page 26 for the quantity O_{K'K}^I
    """
    Olap_I = zeros((len(Km),len(Km)), dtype=float)
    for i in range(len(Km)):
        Olap_I[i,i] = 1 - 4*pi*RMuffinTin**3/(3.*Vol)             # if K=K' => delta(K-K')-4*pi*S^3/V_cell*(lim_dK->0 j1(dKS)/(dKS))
        for j in range(i+1, len(Km)):
            KKl = sqrt(dot(Km[i]-Km[j],Km[i]-Km[j]))              # |K'-K|
            fbessel = special.sph_jn(1,KKl*RMuffinTin)[0][1]      # j1(|K'-K|S)
            Olap_I[i,j] = -4*pi*RMuffinTin**2*fbessel/(KKl*Vol)   # -4*pi*S**2*j1(|K'-K|S)/(V_cell*|K'-K|)
            Olap_I[j,i] = Olap_I[i,j]                             # symmetrize

    return Olap_I

def Wave(Z, Enu, R0, Veff):
    """Solves the SCH Eq for Psi(Enu) and its energy derivative
       Returns logarithmic derivative, Psi(l,E) and its energy derivative
       
       Please  see Eq.38-39 on page 22 for definition, and Eq.49 on page 26
        for S*Psi'(S)/Psi(S) and S*dPsi'(S)/dPsi(S)
    """
    def startSol(Z, l, r):
        "good choice for starting Numerov algorithm"
        return r**(l+1)*(1-Z*r/(l+1))

    logDer=[]
    Psi_l=[]
    Psip_l=[]
    for l in range(len(Enu)):

        # Computes Psi=u/r
        crhs = CRHS(Enu[l], l, R0, Veff)
        crhs[0]=0
        ur = Numerov(crhs, (R0[-1]-R0[0])/(len(R0)-1.), 0.0, startSol(Z,l,R0[1]))
        
        ur *= 1/sqrt( integrate.simps(ur*ur, R0) )  # normalization
        Psi_l.append( ur/R0 ) # storing Psi
        Psi_l[-1][0] = extrapolate(R0[0], R0[1], R0[2], ur[1]/R0[1], ur[2]/R0[2])
        
        # For energy derivative of Psi' = urp/r
        inhom = -2*ur
        urp = NumerovGen(crhs, inhom, (R0[-1]-R0[0])/(len(R0)-1.), 0.0, startSol(Z,l,R0[1]))

        # Energy derivative should be orthogonal
        alpha = integrate.simps(ur*urp, R0)
        urp -= alpha*ur   # can subtract any solution of the homogenious part of the equation, which is ur!
        Psip_l.append( urp/R0 ) # storing Psip'
        Psip_l[-1][0] = extrapolate(R0[0], R0[1], R0[2], urp[1]/R0[1], urp[2]/R0[2])
        
        # <\cdot{\psi}|\cdot{\psi}>
        PsipPsip = integrate.simps(urp*urp, R0)
        
        # Computes the logarithmic derivative
        v1 = crhs[-1]*ur[-1]
        v0 = crhs[-2]*ur[-2]
        w1 = crhs[-1]*urp[-1]+inhom[-1]
        w0 = crhs[-2]*urp[-2]+inhom[-2]
        dh = R0[2]-R0[1]
        dudr  = (ur[-1]-ur[-2])/dh + 0.125*dh*(3*v1+v0)
        dupdr = (urp[-1]-urp[-2])/dh + 0.125*dh*(3*w1+w0)
        
        dlogPsi = RMuffinTin*dudr/ur[-1] - 1     # R*(dpsi/dr)/psi = (du/dr)/u -1
        dlogPsip = RMuffinTin*dupdr/urp[-1] - 1  # ...
        Psi = ur[-1]/RMuffinTin
        Psip = urp[-1]/RMuffinTin
        
        logDer.append( (Psi, Psip, dlogPsi, dlogPsip, PsipPsip) )
        
    return (logDer, Psi_l, Psip_l)

def FindCoreStates(core, R, Veff, Z, fraction=4.):
    "Finds all core states"
    def root(Ex, l, R, Veff):
        "For searching the core bound states"
        rhs = CRHS(Ex, l, R, Veff)
        h = (R[-1]-R[0])/(len(R)-1.)
        u = Numerov(rhs, h, R[0]*exp(-R[0]), R[1]*exp(-R[1]))
        extraplt = u[-2]*(2+h**2*rhs[-2])-u[-3]
        return u[-1]

    coreRho = zeros(len(R), dtype=float)
    coreE = 0
    coreZ = 0

    states=[]
    for l in range(len(core)):
        n=0                           # number of states found
        E = -0.5*Z*Z/(l+1)**2-3.      # here we starts to look for zero
        dE = abs(E)/fraction          # the length of the first step 
        decrease = abs(E)/(abs(E)-dE) # the step will decrease to zero. Check the formula!
        v0 = root(E, l, R, Veff)      # starting value
        while E<0 and n<core[l]:      # we need ncore[l] bound states
            E += dE
            v1 = root(E, l, R, Veff)
            if v1*v0<0:
                Energy = optimize.brentq(root, E-dE, E, args=(l, R, Veff)) # find bound state by root-searching
                # Density
                rhs = CRHS(Energy, l, R, Veff) # now we will compute the radial wave for this state
                u = Numerov(rhs, (R[-1]-R[0])/(len(R)-1.), R[0]*exp(-R[0]), R[1]*exp(-R[1])) # Radial wave
                drho = u*u  # radial charge
                norm = abs(integrate.simps(drho, R ))
                drho *= 1./(norm*4*pi*R**2)
                
                coreRho += drho * (2*(2*l+1.))
                coreE   += Energy*(2*(2*l+1.))
                coreZ   += 2*(2*l+1)
                states.append( (n,l,Energy) )
                n += 1
            dE/=decrease
            v0 = v1

    print '   Found core states for (n,l)=[',
    for state in states:
        print '(%d,%d)' % state[:2],
    print '] E=[',
    for state in states:
        print '%f,' % state[2],
    print ']'
    
    return (coreRho[::-1], coreE, coreZ, states)

def ComputeEigensystem(k, Km, Olap_I, Enu, logDer, RMuffinTin, Vol, VKSi=0):
    """The main part of LAPW algorithm: Implements valence H[K,K'] and O[K,K'] and diagonalizes them.
       Implements all equations on page 26 and page 30.
       The output are energy bands, eigenvectors and weight functions which can be used to compute
       electronic charge in real space.
    """
    def dlog_bessel_j(lmax, x):
        """Calculates logarithmic derivative of the spherical bessel functions
           It returns three quantities:
             (x*d/dx log(j_l(x)),  j_l(x), the product of the first two)
           for l up to lmax
           The last entry is important for singular cases: when x is zero of bessel function. In this case
           the logarithmic derivative is diverging while j*dlog(j(x))/dx is not
        """
        if (fabs(x)<1e-5):
            return [(l, x**l/misc.factorial2(2*l+1), l*x**l/misc.factorial2(2*l+1)) for l in range(lmax+1)]
        else:
            (jls, djls) = special.sph_jn(lmax,x) # all jl's and derivatives for l=[0,...lmax]
            return [(x*djls[l]/jls[l], jls[l], x*djls[l]) for l in range(lmax+1)]


    # Here we prepare coefficients a and b defined in Eq. 60-61 on page 29.
    #   a[K,l] = dot{psi} d/dr j_l - j_l d/dr dot{psi}
    #   b[K,l] = j_l d/dr psi - psi d/dr j_l
    #   PP    -- <\cdot{psi}|\cdot{psi}>
    #   see Eq.47,48 on page 25 for definition.
    a_lk = zeros((len(Km), len(Enu)), dtype=float)
    b_lk = zeros((len(Km), len(Enu)), dtype=float)
    PP = array([logDer[l][4] for l in range(len(Enu))])
    for iK,K in enumerate(Km):
        Dl_jl = dlog_bessel_j(len(Enu)-1, sqrt(dot(k+K,k+K))*RMuffinTin)
        for l in range(len(Enu)):
            (Psi, Psip, dlogPsi, dlogPsip, PsipPsip) = logDer[l]
            (Dl, jl, jlDl) = Dl_jl[l]
            # Dl = x/j_l * (d j_l/dx) ; Dl_jl = x * d j_l/dx
            a_lk[iK,l] =  jl * Psip/RMuffinTin * (Dl - dlogPsip)
            b_lk[iK,l] = -jl * Psi/RMuffinTin  * (Dl -  dlogPsi)
            
    # This part of code is too slow in Python, hence was recoded in C++
    # It computes the argument of the Legendre polynomial (see Eq.59 on page 29).
    #     argums(K,K') = (k+K)*(k+K')/(|k+K| * |k+K'|)
    #     qv(K)  = k+K
    #     qvs(K) = |k+K|
    # where K is reciprocal vector and k is momentum vector.
    codeArg="""
         for (int iK=0; iK<Km.extent(0); iK++){ // computing qvs <= |k+K|
            for (int i=0; i<3; i++) qv(iK,i) = Km(iK,i)+k(i);
            qvs(iK) = sqrt( qv(iK,0)*qv(iK,0) + qv(iK,1)*qv(iK,1) + qv(iK,2)*qv(iK,2) );
         }
         for (int iK=0; iK<Km.extent(0); iK++){
            for (int jK=0; jK<Km.extent(0); jK++){
               double qvqv=0;
               for (int i=0; i<3; i++) qvqv += qv(iK,i)*qv(jK,i);
               if (qvs(iK)*qvs(jK)==0) argums(iK,jK)=1.;
               else argums(iK,jK) = qvqv/(qvs(iK)*qvs(jK));
            }
         }
    """
    qv = zeros((len(Km),3), dtype=float)
    qvs = zeros(len(Km), dtype=float)
    argums = zeros((len(Km), len(Km)), dtype=float)  # becomes (k+K).(k+K')/(|k+K||k+K'|)
    weave.inline(codeArg, ['Km', 'k', 'qv', 'qvs', 'argums'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    
    # This part of code is too slow in Python even when using special.Legendre
    # It computes the Legendre polynomial for all values of argums(K,K') precomputes above
    # and for all l up to lmax=len(Enu)
    # The few lowest order Legendre polynomials are precomputes and recursion is used only for high order l.
    #   Computes:  Leg(K,K',l) = P_l(argums(K,K'))
    codeLegendre="""
      for (int iK=0; iK<argums.shape()[0]; iK++){
         for (int jK=0; jK<argums.shape()[1]; jK++){
            double x=argums(iK,jK);
            double x2 = x*x;
            Leg(iK,jK,0) = 1;
            if (lmax>=1)  Leg(iK,jK,1) = x;
            if (lmax>=2)  Leg(iK,jK,2) = 1.5*x2-0.5;
            if (lmax>=3)  Leg(iK,jK,3) = x*(2.5*x2-1.5);
            if (lmax>=4)  Leg(iK,jK,4) = 0.375*(1-10*x2*(1-1.1666666666666667*x2));
            if (lmax>=5)  Leg(iK,jK,5) = 1.875*x*(1-4.66666666666666667*x2*(1-0.9*x2));

            for (int l=6; l<=lmax; l++){
                double p0 = 0.375*(1-10*x2*(1-1.1666666666666667*x2));
                double p1 = 1.875*x*(1-4.66666666666666667*x2*(1-0.9*x2)); 
                double p2=0;
                for (int i=6; i<=l; i++){
                  p2 = ((2*i-1)*x*p1-(i-1)*p0)/i;
                  p0=p1;
                  p1=p2;
                }
                Leg(iK,jK,l) = p2;
            }
         }
      }
    """
    lmax = len(Enu)-1
    Leg = zeros((len(Km),len(Km),len(Enu)), dtype=float)
    weave.inline(codeLegendre, ['argums', 'lmax', 'Leg'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    
    # This part of code is too slow in Python
    # Implements Eq. 55-59 on page 29
    # It computes the Hamiltonian and Overlap in Muffin-Thin and interstitials
    # All necessary arrays were precomputes above.
    codeHO="""
        for (int iK=0; iK<Ham.extent(0); iK++){
            for (int jK=0; jK<Ham.extent(1); jK++){
                double olapMT=0;
                double hamMT=0;
                for (int l=0; l<Enu.size(); l++){
                    double Plk = C0n(l)*Leg(iK,jK,l);
                    double a_a = a_lk(iK,l) * a_lk(jK,l);
                    double a_b = a_lk(iK,l) * b_lk(jK,l);
                    double b_a = b_lk(iK,l) * a_lk(jK,l);
                    double b_b = b_lk(iK,l) * b_lk(jK,l);
    
                    double olap = a_a + b_b*PP(l);
                    olapMT += Plk * olap;
                    hamMT  += Plk * ( 0.5*(b_a + a_b) + olap*Enu(l) );

                    WK0(l,iK,jK)  = Plk *  a_a;
                    WK1(l,iK,jK) = Plk * (b_a + a_b);
                    WK2(l,iK,jK) = Plk *  b_b;
                }
                Olap(iK,jK) = olapMT + Olap_I(iK,jK);
                Ham(iK,jK) = ( 0.25*(qvs(iK)*qvs(iK) + qvs(jK)*qvs(jK)) + VKSi )*Olap_I(iK,jK) + hamMT;
            }
        }
    """
    Olap = zeros((len(Km), len(Km)), dtype=float)
    Ham  = zeros((len(Km), len(Km)), dtype=float)
    WK0 = zeros((len(Enu), len(Km), len(Km)), dtype=float)
    WK1 = zeros((len(Enu), len(Km), len(Km)), dtype=float)
    WK2 = zeros((len(Enu), len(Km), len(Km)), dtype=float)
    Enu = array(Enu)
    C0n = (2*arange(len(Enu)) + 1) * pi * RMuffinTin**4/Vol
    weave.inline(codeHO, ['Olap', 'Ham', 'WK0', 'WK1', 'WK2', 'Enu', 'Leg', 'qvs', 'PP', 'Olap_I', 'VKSi', 'a_lk', 'b_lk', 'C0n'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    # Diagonalizes the LAPW Hamiltonian
    #Ek, Ar = symeig.symeig(Ham, Olap, type=1) # symmetric generalized eigenvalue problem
    Ek, Ar = linalg.eigh(Ham, Olap, type=1) # symmetric generalized eigenvalue problem

    
    #print matrix(Ar).T * matrix(Ham) * matrix(Ar)
    
    # Calculation of weights for valence density
    # Implements the weight functions Eqs.70-73 on page 33.
    Ar = matrix(Ar)
    w0 = zeros((len(Enu),len(Ar)),dtype=float)
    w1 = zeros((len(Enu),len(Ar)),dtype=float)
    w2 = zeros((len(Enu),len(Ar)),dtype=float)
    for l in range(len(Enu)):
        tw0 = Ar.H * matrix(WK0[l]) * Ar
        tw1 = Ar.H * matrix(WK1[l]) * Ar
        tw2 = Ar.H * matrix(WK2[l]) * Ar
        w0[l,:] = array([tw0[p,p] for p in range(len(Ar))])
        w1[l,:] = array([tw1[p,p] for p in range(len(Ar))])
        w2[l,:] = array([tw2[p,p] for p in range(len(Ar))])
        #w0[l,:] = trace(Ar.H * matrix(WK0[l]) * Ar)
        #w1[l,:] = trace(Ar.H * matrix(WK1[l]) * Ar)
        #w2[l,:] = trace(Ar.H * matrix(WK2[l]) * Ar)
        
    twi = Ar.T * matrix(Olap_I) * Ar  # Here we approximate the charge in the interstitals by a constant, which gives the weight equal to Overlap^I.
    wi = array([twi[p,p] for p in range(len(Ar))])

    return (Ek, Ar, w0, w1, w2, wi)



def rootChemicalPotential(mu, Ek, wkp, Zval, beta=50.):
    " Computes valence density to find root for the chemical potential "
    code_mu="""
    double Zt=0;
    for (int ik=0; ik<Ek.shape()[0]; ik++){
        for (int p=0; p<Ek.shape()[1]; p++){
            double x = beta*(Ek(ik,p)-mu);
            double ferm = abs(x) < 100 ? 1/(exp(x)+1) : (x<0 ? 1 : 0);
            Zt += wkp(ik) * ferm;
        }
    }
    return_val = Zt;
    """
    Zt = weave.inline(code_mu, ['mu', 'Ek', 'beta', 'wkp'],type_converters=weave.converters.blitz, compiler = 'gcc')
    return 2*Zt-Zval


def ComputeMTDensity(mu, Ek, wkp, w0, w1, w2, Psi_l, Psip_l, beta=50.):
    """Given the coefficients Eqs.70-73 on page 33, it computes the valence charge
       given the chemical potential mu. The Eq. 75 on page 34 descibes the algorithm.
       Note that the radial solution of the SCH equation is necessary together with its energy derivative
    """
    code_sum="""
    double dw0, dw1, dw2; dw0=dw1=dw2=0;
    for (int p=0; p<ek.shape()[0]; p++){
        double x = beta*(ek(p)-mu);
        double ferm;
        ferm = abs(x) < 100 ? 1/(exp(x)+1) : (x<0 ? 1 : 0);
        dw0 += w0k(p)*ferm;
        dw1 += w1k(p)*ferm;
        dw2 += w2k(p)*ferm;
    }
    //return_val = make_tuple( dw0, dw1, dw2 ); // does not yet work
    py::tuple results(3);
    results[0] = dw0; results[1] = dw1; results[2] = dw2;
    return_val = results;
    """
    nlmax = len(w0[0])
    wgh=zeros((nlmax,3), dtype=float)
    for l in range(nlmax):
        for ik in range(shape(Ek)[0]):
            w0k = array(w0[ik][l,:]) # weight for all bands
            w1k = array(w1[ik][l,:])
            w2k = array(w2[ik][l,:])
            ek = array(Ek[ik])      # eigenvalues for all bands
            dws = weave.inline(code_sum, ['mu', 'ek', 'beta', 'w0k', 'w1k', 'w2k'],type_converters=weave.converters.blitz, compiler = 'gcc')
            wgh[l,:] += array(dws) * wkp[ik] # accumulates weight for each l and [psi*psi,dot{psi}*psi,dot{psi}*dot{psi}]
        #print "%3d %20.10f %20.10f %20.10f" % (l, wgh[l,0], wgh[l,1], wgh[l,2])

    # Implements Eq.75 on page 34.
    nR = len(Psi_l[0])
    MTRho = zeros(nR, dtype=float)
    for l in range(nlmax):
        for ir in range(nR):
            MTRho[ir] += wgh[l,0]*Psi_l[l][ir]**2 + wgh[l,1]*Psi_l[l][ir]*Psip_l[l][ir] + wgh[l,2]*Psip_l[l][ir]**2
            
    MTRho *= 2/(4*pi)  # 2 due to spin

    return MTRho

def ComputeInterstitialCharge(mu, Ek, wi, wkp, beta=50.):
    " Interstitial charge Eq.64 page 31"
    def ferm(x):
        if x>100: return 0.0
        if x<-100: return 1.0
        return 1/(exp(x)+1.)

    sIntRho=0  #  total interstitial charge
    for ik in range(len(Ek)):
        dsum = 0.0
        for p in range(len(Ek[ik])):
            dsum += ferm( (Ek[ik,p]-mu)*beta )*wi[ik][p]
        sIntRho += dsum*wkp[ik]
    sIntRho *= 2 # due to spin
    return sIntRho


#######################################################
###     Atomic charge for a good starting guess   #####
#######################################################
def Atom_cmpb(x,y):
    "Comparison function for sorting of bound states"
    if abs(x[2]-y[2])>1e-4:
        return cmp(x[2],y[2])
    else:
        return cmp(x[1],y[1])

def Atom_ChargeDensity(states, R, Veff, Z):
    "Computes electron charge density, given the bound states and Z"
    rho = zeros(len(R), dtype=float)
    N = 0
    Ebs = 0
    for state in states:
        l = state[1]
        E = state[2]
        rhs = CRHS(E, l, R, Veff)
        h = (R[-1]-R[0])/(len(R)-1.)
        u = Numerov(rhs, h, R[0]*exp(-R[0]), R[1]*exp(-R[1]))
        u /= sqrt(abs(sum(u)*h))   # To avoid overflow
        
        u2 = u**2
        norm = abs(integrate.simps(u2, R))
        u2 *= 1./norm

        dN = 2*(2*l+1.)
        if N+dN<Z:
            ferm = 1.
        else:
            ferm = (Z-N)/dN
        drho = u2[:]*dN*ferm/(4*pi*R**2)
        
        rho += drho

        N += dN
        Ebs += E * ferm*dN
        if N>=Z: break

    rho[-1] = extrapolate(0.0, R[-2], R[-3], rho[-2], rho[-3])
    return (rho, Ebs)

def Atom_charge(Z, core, mix=0.3, RmaxAtom=10., Natom=3001, precision=1e-5, Nitt=100):
   """ Computes Atomic electronic density and atomic Energy
   Input:
      Z             --  Nucleolus charge
      core          --  States treated as core in LAPW (example: [3,2,0]  # 1s,2s,3s, 1p,2p, no-d)
      mix           --  Mixing parameter for density
      RmaxAtom      --  The end of the radial mesh (maximum r)
      Natom         --  Number of points in radial mesh
      precision     --  How precise total energy we need      
      Nitt          --  Maximum number of itterations
      """

   XC = excor.ExchangeCorrelation(3)  # Exchange correlations class; WVN seems to be the best (look http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html)

   R0 = linspace(1e-10, RmaxAtom, Natom) # Radial mesh
   Ra = R0[::-1]                         # Inverse radial mesh. Note that we start shooting from the MT-boundary down to r=0 for numeric stability.

   Veff = -ones(len(Ra), dtype=float)/Ra # Hydrogen-atom type nucleous

   catm = [c + 1 for c in core]          # We add one more state to core to get atomic states
   
   Etot_old = 0
   # Finds bound states
   (coreRho, coreE, coreZ, states) = FindCoreStates(catm, Ra, Veff, Z)
   # Sorts them according to energy
   states.sort(Atom_cmpb)
   # Computes charge
   (rho, Ebs) = Atom_ChargeDensity(states, Ra, Veff, Z)
   rho=rho[::-1]

   for itt in range(Nitt):

       # Here we have increasing R ->

       # Hartree potential
       UHartree = SolvePoisson(Z, R0, rho)
       # Adding exchange-correlation part
       Vxc = [XC.Vx(rsi)+XC.Vc(rsi) for rsi in rs(rho)]
       ExcVxc = [XC.EcVc(rsi) + XC.ExVx(rsi) for rsi in rs(rho)]
       Veff = (UHartree - Z)/R0 + Vxc
       Veff=Veff[::-1]  # Turn around for later calculations.

       # Here we have decreasing R <-
       
       # Finds bound states
       (coreRho, coreE, coreZ, states) = FindCoreStates(catm, Ra, Veff, Z)
       # Sorts them according to energy
       states.sort(Atom_cmpb)
       # Computes charge
       (nrho, Ebs) = Atom_ChargeDensity(states, Ra, Veff, Z)

       # Total energy
       pot = (ExcVxc*R0**2-0.5*UHartree*R0)*nrho[::-1]*4*pi

       Etot = integrate.simps(pot, R0) + Ebs
       Ediff = abs(Etot-Etot_old)
       
       print '   %d) Etot=%f Eband=%f Ediff=%f' % (itt, Etot, Ebs, Ediff)
       
       # Mixing
       rho = mix*nrho[::-1] + (1-mix)*rho
       Etot_old = Etot

       if Ediff < precision: break

   return (R0, rho)
   

#################################
### Small utility functions   ###
#################################
def extrapolate(x, x0, x1, f0, f1):
    "linear extrapolation"
    return f0 + (f1-f0)*(x-x0)/(x1-x0)

def rs(rh):
    "rs from density -> an electron radius that corresponds to density"
    return (3/(4*pi*rh))**(1/3.)

DEFAULT_COLOR = '\033[0m'
RED = '\033[31;1m'
GREEN = '\033[32;1m'
BLUE = '\033[34;1m'
YELLOW = '\033[33;1m'

if __name__ == '__main__':
    ###################################
    # Start input parameters
    Z=29                     # Number of electrons in the atom
    LatConst = 6.8219117     # Lattic constant
    nkp = 6                  # Number of k-points in 1BZ: (nkp x nkp x nkp)
    #### Core states ##############
    core = [3,2,0]  # [#s-states, #p-states, #d-states], i.e., 1s,2s,3s, 1p,2p, no-d
    #### Linearization energies ###
    # Most of high energy partial waves should be centered around mu
    # In general, it is a good idea to determine them self-consistently
    Enu = [0.11682, 0.18794, 0.211145, 0.3, 0.3, 0.3]
    N = 1001                 # Number of points in radial mesh
    beta=50.                 # To compute chemical potential we take finite inverse temperature
    mu_mm = [0.0, 1.0]       # Chemical potential is searched between mu_mm[0] and mu_mm[1]
    CutOffK=3.5              # Largest lengt of reciprocal vectors K (only shorter vec. are taken into account)
    DRho = 1e-3              # Convergence criteria for electronic density Rho
    Nitt = 100               # Maximum number of itterations
    mixRho = 0.3             # Linear mixing parameter for charge
    Nkplot = 200             # Number of k-points for plotting bands
    plotMM = [-1.,0.1]       # Bands will be plotted in this energy range
    # End inout parameters 
    ########################
    
    # Core number of electrons
    Zcor = sum([2*(2*l+1)*nl for l,nl in enumerate(core)])
    # Valence number of electrons
    Zval = Z-Zcor
    
    print "Z core=", Zcor, " and Zval=", Zval
    
    # Creates atomic charge to have a good starting point
    (Atom_R0, Atom_rho) = Atom_charge(Z, core, 0.3)
    AtomRhoSpline = interpolate.UnivariateSpline(Atom_R0, Atom_rho, s=0)
    
    # Exchange correlations class; WVN seems to be the best (look http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html)    
    XC = excor.ExchangeCorrelation(3)
    
    # Generates and stores momentum points
    fcc = FccLattice(LatConst)                  # Information about lattice
    RMuffinTin = fcc.RMuffinTin()               # Muffin-Tin radius choosen such that spheres touch
    VMT = 4*pi*RMuffinTin**3/3.                 # Volume of MT
    Vinter = fcc.Volume-VMT                     # Volume of the interstitial region
    print "Muffin-Tin radius =", RMuffinTin
    print "Volume of the MT sphere    =", VMT
    print "Volume of the unit cell    =", fcc.Volume
    print "Volume of the interstitial =", Vinter
    fcc.GenerateReciprocalVectors(4, CutOffK)   # Reciprocal bravais lattice is builded, K points taken into account only for |K|<CutOff
    fcc.ChoosePointsInFBZ(nkp,0)                # Chooses the path in the 1BZ or the k-points in the irreducible 1BZ

    # Radial mesh --  only linear mesh can be used in connection to Numerov algorithm.
    R0 = linspace(0, RMuffinTin, N)
    R0[0]=1e-10
    R = R0[::-1] # mesh from RMT down to zero.
    
    # Interstital overlap does not change through iterations
    Olap_I = ComputeInterstitialOverlap(fcc.Km, RMuffinTin, fcc.Volume)

    # We interpolate atomic charge on the new mesh within Muffin-Tin sphere
    TotRho = AtomRhoSpline(R0)
    
    for itt in range(Nitt):  # self-consistent loop
        
        print '%d) Preparing potential' % itt
        UHartree = SolvePoisson(Z, R0, TotRho)
        # Adding exchange-correlation part
        Vxc = [XC.Vx(rsi)+XC.Vc(rsi) for rsi in rs(TotRho)]
        
        nVeff = (UHartree - Z)/R0 + Vxc
        zeroMT = nVeff[-1]  # New MT zero
        nVeff -= zeroMT
        print '   Muffin-Tin zero is ', zeroMT

        Veff = nVeff
        
        (logDer, Psi_l, Psip_l) = Wave(Z, Enu, R0, Veff)
        
        (coreRho, coreE, coreZ, core_states) = FindCoreStates(core, R0[::-1], Veff[::-1], Z)
        print '   coreZ=', coreZ, 'coreE=', coreE
        
        # This is the main loop over all k-points
        Ek=[]; w0=[]; w1=[]; w2=[]; wi=[]
        for ik,k in enumerate(fcc.kp):
            (tEk, tAr, tw0, tw1, tw2, twi) = ComputeEigensystem(k, fcc.Km, Olap_I, Enu, logDer, RMuffinTin, fcc.Volume)
            Ek.append(tEk); w0.append(tw0); w1.append(tw1); w2.append(tw2); wi.append(twi);
        Ek = array(Ek)

        # New chemical potential
        mu = optimize.brentq(rootChemicalPotential, mu_mm[0], mu_mm[1], args=(Ek, fcc.wkp, Zval, beta))
        print GREEN, 'New chemical potential is', mu, DEFAULT_COLOR
        
        MTRho = ComputeMTDensity(mu, Ek, fcc.wkp, w0, w1, w2, Psi_l, Psip_l, beta)
        nTotRho = MTRho + coreRho
        
        sMTRho = integrate.simps(MTRho*R0**2*(4*pi), R0)
        sIntRho = ComputeInterstitialCharge(mu, Ek, wi, fcc.wkp, beta)
        sCoreRho = integrate.simps(coreRho*R0**2*(4*pi), R0)
        
        print '   Zval=', Zval, '~', sMTRho+sIntRho
        print '   Weght in the MT sphere =', sMTRho, 'and in the interstitials =', sIntRho, 'and in core =', sCoreRho
        
        renorm = Z/(sMTRho+sIntRho+sCoreRho)
        print '   Total charge found=', sMTRho+sIntRho+sCoreRho, 'should be', Z, '-> renormalizing by', renorm
        nTotRho *= renorm

        DiffRho = integrate.simps(abs(nTotRho-TotRho), R0)
        print BLUE, 'Electron density difference=', DiffRho, DEFAULT_COLOR 
        if (DiffRho<DRho): break
        
        TotRho = mixRho * nTotRho + (1-mixRho)*TotRho
    

    # Plotting bands
    fcc.ChoosePointsInFBZ(Nkplot, type=1) 
    
    Ek=[]
    for ik,k in enumerate(fcc.kp):
        (tEk, tAr, tw0, tw1, tw2, twi) = ComputeEigensystem(k, fcc.Km, Olap_I, Enu, logDer, RMuffinTin, fcc.Volume)
        Ek.append(tEk)
    Ek = array(Ek)


    for i in range(shape(Ek)[1]):
        if max(Ek[:,i])-mu > plotMM[0] and min(Ek[:,i])-mu < plotMM[1]:
            plot(Ek[:,i]-mu, 'k-', lw=2)

    plot([0,len(Ek)],[0,0], 'k:')  # chemical potential line
    ax=axis()

    xlabs = [p[1] for p in fcc.Points]
    labs  = [p[0] for p in fcc.Points]
    xticks(xlabs, labs)
    
    for ix,x in enumerate(xlabs):
        plot([x,x], [ax[2],ax[3]], 'k:')
        
    axis([xlabs[0], xlabs[-1], ax[2], ax[3]])
    show()

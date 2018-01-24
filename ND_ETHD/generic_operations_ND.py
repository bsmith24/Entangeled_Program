#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import math
import methods_ND

#"""
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
#"""

def harmonic_oscillator(q, w, m, dim):

    """
    Defines the potential for the Harmonic Oscillator
    This potential can be expressed in N dimensions
 
    1D: V(q[0]) = 0.5*m*w*w*(q[0]*q[0])
    2D: V(q1,q2) = 0.5*m*w*w(q1*q1 + q2*q2)

    Params in:  q = list of trajectory positions
                w = potential specfifc constant 
                m = list of masses 
                dim - # of dof for each trajectory
    """    

    # Initilize list for potential energy and force
    # 1D
    # pot will be [0.0, 0.0], where pot[0] holds the sum of classical contributions
    # f will be a list: [0,1,2,...,N] of forces, where each element of the list is the force for a particular particle

    # 2D
    # pot will be [ [0.0, 0.0] , [0.0, 0.0]  ] , where pot[0] holds the sum classical contributions
    # the form is: pot = [ [1st_dof_classic is here , 2nd_dof_classic is here] , [0.0, 0.0] ]
    # f will be a list of two component lists: [ [0,0],[1,1],[2,2],...,[N,N] ] 
    
    N = len(q)
    pot = [ [0.0]*dim , [0.0]*dim ]
    f = []

    """
    print "pot = ", pot
    print "f = ", f
    print "\n"
    sys.exit(0)
    """
    for i in xrange(N):
        f.append( [0.0]*dim )    
        for j in xrange(dim):

            pot[0][j] = pot[0][j] +  0.5*m[i]*w*w*q[i][j]*q[i][j]

            # Recall q[particle_index][dimension] -> q[0][1] = first particle, 2nd dof
            f[i][j] = -m[i]*w*w*q[i][j] 

    """
    print "pot = ", pot
    print "f = ", f
    print "\n"
    sys.exit(0)
    """

    return pot, f

def double_well_potential(q, dim):

    """
    This potential can be expressed in N dimensions
    1D expression: V(q) = A * ( 0.25*q^4 - 0.5*q^2 )

    Params in:  q = list of trajectory positions
                dim - # of dof for each trajectory
    """

    C = 1.0
    a, b = 0.25, 0.5

    N = len(q)
    pot = [ [0.0]*dim , [0.0]*dim ]
    f = []

    q2, q3, q4 = [], [], []
    for i in xrange(N):

        f.append( [0.0]*dim )
        q2.append( [0.0]*dim )
        q3.append( [0.0]*dim )
        q4.append( [0.0]*dim )    

        for j in xrange(dim):

            q2[i][j] = q[i][j]*q[i][j]
            q3[i][j] = q2[i][j]*q[i][j]
            q4[i][j] = q3[i][j]*q[i][j]

            pot[0][j] = pot[0][j] + C*( a*q4[i][j] - b*q2[i][j] )

            # Recall q[particle_index][dimension] -> q[0][1] = first particle, 2nd dof
            f[i][j] =  -C*( 4.0*a*q3[i][j] - 2.0*b*q[i][j] )

    """
    print "pot = ", pot
    print "f = ", f
    print "\n"
    sys.exit(0)
    """

    return pot, f


def Martens_Type1(q, dim):

    """
    This potential is coded strictly for two dimensions:

    V(q1,q2) = Va*sech^2(2*q1) + 0.5*Vb*q2*q2
    sech(q) = 1/(math.cosh(q)
    sech^2(q) = (1/math.cosh(q))*(1/math.cosh(q))

    Params in:  q = list of trajectory positions
                dim - # of dof for each trajectory 9stricktly 2 for this potential)
    """

    if dim != 2:
        print "\n", "ERROR"
        print "Must use 2 degrees of freedom in order to use model potential: Marten's Type1"
        print "Please check model potential type or try adjusting the variable dim to 2"
        sys.exit(0)
  
    else:

        # Define Constants
        Va, Vb = 0.00625, 0.0106

        N = len(q)
        pot = [ [0.0]*dim , [0.0]*dim ]
        f = []
        for i in xrange(N):

            f.append( [0.0]*dim )    

            sech =  1.0/math.cosh(2.0*q[i][0])
            sech2 = sech*sech

            pot[0][0] = pot[0][0] + Va*sech2
            pot[0][1] = pot[0][1] + 0.5*Vb*q[i][1]*q[i][1] 

     
            # Recall q[particle_index][dimension] -> q[0][1] = first particle, 2nd dof
            f[i][0] = 0.025*math.tanh(2.0*q[i][0])*sech2
            f[i][1] = - Vb*q[i][1]


        """
        print "pot = ", pot
        print "f = ", f
        print "\n"
        sys.exit(0)
        """

        return pot, f


def Martens_Type2(q, dim):

    """
    Currently not working - 1/23/17 
    Declared working on - 

    This potential is strictly coded for two dimensions:

    V(q1,q2) = Va*sech^2(2*q1) + 0.5*Vb*(q2 + Vc*(q1^2 - 1.0))^2
    sech(q) = 1/(math.cosh(q)
    sech^2(q) = (1/math.cosh(q))*(1/math.cosh(q))

    Params in:  q = list of trajectory positions
                dim - # of dof for each trajectory
    """

    if dim != 2:
        print "\n", "ERROR"
        print "Must use 2 degrees of freedom in order to use model potential: Marten's Type2"
        print "Please check model potential type or try adjusting the variable dim to 2"
        sys.exit(0)
  
    else:

        # Define Constants
        Va, Vb, Vc = 0.00625, 0.0106, 0.4

        N = len(q)
        pot = [ [0.0]*dim , [0.0]*dim ]
        f = []
        for i in xrange(N):

            f.append( [0.0]*dim )    

            x = q[i][0]
            y = q[i][1]

            sech =  1.0/math.cosh(2.0*x)
            sech2 = sech*sech

            # Here, pot[0][0] will carry the non coupled part of the potential
            # pot[0][1] will carry the coupled part of the potential

            pot[0][0] = pot[0][0] + Va*sech2
            pot[0][1] = pot[0][1] + 0.5*Vb*(y + Vc*(x*x - 1.0))**2.0 

            # Recall q[particle_index][dimension] -> q[0][1] = first particle, 2nd dof
            f[i][0] = - (0.00848*x*(0.4*(x*x - 1.0) + y) - 0.025*math.tanh(2.0*x)*sech2)
            f[i][1] = - (0.00424*x*x + 0.0106*y - 0.00424)

        """
        print "pot = ", pot
        print "f = ", f
        print "\n"
        sys.exit(0)
        """

        return pot, f

def compute_kin(p,m, dim):

    """
    Defines the kinetic energy for a Hamiltonian
    1D: T(p[0]) = 0.5*(p[0]*p[0])/m
    2D: T(p[i][0],p[i][1]) = 0.5*(p[i][0]*p[i][0] + p[i][1]*p[i][1])/m

    Params in:
    p = momentum of particle
    m = list of masses 
    """    

    # The kinetic energy is structured like the potential energy for 2D, except we will not have the additional list for quantum forces
    # The kinetic energy is structured as follows: kin = [ [0.0,0.0] ]
    # In the form [ [1st_dof_kinetic_here , 2nd_dof_kinetic_here] ]

    N = len(p)
    kin = [ [0.0]*dim ]

    for i in xrange(N):
        for j in xrange(dim):

            # Kinetic: for the x dimension of each particle
            kin[0][j] = kin[0][j] + 0.5*p[i][j]*p[i][j]/m[i]

    return sum(kin[0])

def propagate_p(p, f, dim, dt):
# Returns the updated momenta

    N = len(p)
    for i in range(N):
        for j in xrange(dim):

            p[i][j] = p[i][j] + f[i][j]*dt

    return p

def propagate_q(q, p, m, dim, dt):
# Returns the updated positions

    N = len(q)
    for i in range(N):
        for j in xrange(dim):
        
            q[i][j] = q[i][j] + p[i][j]*dt/m[i]

    return q

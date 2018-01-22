#=============================================================#
# == Copyright (C) 2017 Brendan A. Smith, Alexey V. Akimov == #
#=============================================================#
import sys
import os
import math
import generic_operations
import methods

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def rescale(q, p, w, m, q_max, approx, model, vel, ent_type):

    # Compute pre-absorption energy:
    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

    T_old = generic_operations.compute_kin(p, m)
    E_old = sum(V)

    q_new, p_new = generic_operations.traj_absorb(q, p, q_max)


    # Compute post-absorption energy:
    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

    T_new = generic_operations.compute_kin(p_new, m)
    E_new = sum(V)

    # Lets compensate for the energy change by rescaling momenta
    # instead of the old kinetic energy (subtract), we have a new one (add)
    # E_old + T_old = E_new + alp^2 * T_new
    # alp^2 = (E_old - E_new + T_old)/T_new 
    
    alp2 = 1.0 
    if T_new>0.0:
        alp2 = (E_old - E_new + T_old)/T_new    # Ask

    alp = 1.0
    if alp2>0.0:

        alp = math.sqrt(alp2)

    if vel == 0:
        alp = 1.0

    N = len(p_new)
    if alp<1.0:
        for i in range(N):        
            p_new[i] =  alp*p_new[i] 


    return q_new, p_new


def propagate_methods(q, p, w, m, f, dt, approx, model, vel, ent_type):
    ############################################################
    # -------- Quantum Propagation: Entangeled Verlet -------- #
    ############################################################

    """

    N = len(q)
    for i in range(N):
        p[i] = p[i] + f[i]*dt*0.5
        q[i] = q[i] + p[i]*dt/m
        V, f = compute_pot_force(q, w, m, approx, model)
        p[i] = p[i] + f[i]*dt*0.5


    """

    generic_operations.propagate_p(p, f, dt*0.5)

    generic_operations.propagate_q(q, p, m, dt)

    q, p = rescale(q, p, w, m, q_max, approx, model, vel, ent_type)

    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

    generic_operations.propagate_p(p, f, dt*0.5)

    return q, p

#############################################################################
def main(q, p, q_grid, p_grid, Nsteps, Nsnaps, m, dt, approx, model, ent_type, vel):

    N = len(q)
    t = 0.0; orig = float(N); barrier = 0.0; QG = len(q_grid); PG = len(p_grid) 
             
    os.system("mkdir p_q_info")
    init_q = open("p_q_info/initial_q.txt","w")
    init_p = open("p_q_info/initial_p.txt","w")
    for i in range(N):
        init_q.write( " %8.5f" % (q[i]) )
        init_q.write("\n")
        init_p.write( " %8.5f" % (p[i]) )
        init_p.write("\n")

    w = generic_operations.compute_omega(q, m, approx, model)
    print "omega = ", w
    
    if model == 1:
        barrier = 0.67
        q_max = 5.0
    if model == 2:
        barrier = 0.0
        q_max = 50.0

    # Initilize forces
    if ent_type == 1:
        V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
    if ent_type == 2:
        V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
    if ent_type == 3:
        V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

    print "Calling RPMG force, see if V[1] Matches"
    print "V[0] = ", V[0]/float(N)
    print "V[1] = ", V[1]/float(N)

    os.system("mkdir energy")
    os.system("mkdir phase_space")
    os.system("mkdir pos_space")
    os.system("mkdir distribution_data")
    os.system("mkdir tunnel")
    os.system("mkdir perturbation")
    os.system("mkdir input")

    e = open("energy/energy.txt", "w")
    r = open("phase_space/phase_space.txt", "w")
    g = open("pos_space/pos_space.txt", "w")
    h = open("distribution_data/dist.txt", "w")
    hh = open("tunnel/tunnel.txt", "w")
    s = open("perturbation/perturbation.txt", "w")

    s.write(" %8.5f  %8.5f" % (t, V[1]/float(N)) ),
    s.write( "\n" )

    V = sum(V)
    T = generic_operations.compute_kin(p, m)
    E = V + T
    
    ### Printing Initial Distribution Information
    for j in range(QG):
        count = 0.00
        for i in range(N):
            if q[i] > q_grid[j-1] and q[i] < q_grid[j]:  
                count += 1.00
            prob = count/float(N)             
        h.write( str(q_grid[j]) + " " + str(prob) + "\n") 

    ###==Printing Stuff==###
    e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T/float(N), V/float(N),E/float(N)) )

    for i in range(N):
        r.write(" %8.5f  %8.5f" % (q[i], p[i]) ),
    r.write( "\n" )

    g.write( " %8.5f" % (t) ),
    for i in range(N):
        g.write( " %8.5f" % (q[i]) )
    g.write( "\n" )

    # Propagate for Nsnaps
    for k in xrange(Nsnaps):

        """
        a = open("input/2D_phase_"+str(k)+".dat","w")
        for j in range(QG):
            for l in range(PG):
                count = 0.00
                for i in range(N):
                    if q[i] > q_grid[j-1] and q[i] <= q_grid[j]:
                       if p[i] > p_grid[l-1] and p[i] <= p_grid[l]:
                           count += 1.00
                prob = count/float(N)
                a.write(" %8.5f  %8.5f  %8.5f" % (q_grid[j], p_grid[l], prob) )
                a.write("\n")
            a.write("\n")
        """

        # Propagate for Nsteps
        for j in range(Nsteps):

#            q, p = propagate_methods(q, p, w, m, f, dt, approx, model, vel, ent_type)

            t = t + dt

#            q, p = rescale(q, p, w, m, q_max, approx, model, vel, ent_type)

            generic_operations.propagate_p(p, f, dt*0.5)
            generic_operations.propagate_q(q, p, m, dt)

            q, p = rescale(q, p, w, m, q_max, approx, model, vel, ent_type)

            if ent_type == 1:
                V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
            elif ent_type == 2:
                V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
            elif ent_type == 3:
                V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)

            generic_operations.propagate_p(p, f, dt*0.5)

            N = len(q)

#            if ent_type == 1:
#                V, f = methods.compute_ETHD_pot_force(q, m, approx, model)
#            if ent_type == 2:
#                V, f = methods.compute_RPMD_pot_force(q, w, m, approx, model)
#            if ent_type == 3:
#                V, f = methods.compute_QHD2_pot_force(q, w, m, approx, model)


        s.write(" %8.5f  %8.5f" % (t, V[1]/float(N)) ),
        s.write( "\n" )

        V = sum(V)
        T = generic_operations.compute_kin(p, m)
        E = V + T

        ###==Printing Stuff==###
        e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, T/float(N), V/float(N),E/float(N)) )

        for i in range(N):
            r.write(" %8.5f  %8.5f" % (q[i], p[i]) ),
        r.write( "\n" )

        g.write( " %8.5f" % (t) ),
        for i in range(N):
            g.write( " %8.5f" % (q[i]) )
        g.write( "\n" )

        count = 0.0
        for i in range(N):
            if q[i] < barrier:
                count = count + 1.0
        count = count/orig
        hh.write( str(t) + " " + str(count) + " " + "\n")

    post_q = open("p_q_info/post_therm_q_dist.txt","w")
    post_p = open("p_q_info/post_therm_p_dist.txt","w")
    for i in range(N):
        post_q.write( " %8.5f" % (q[i]) )
        post_q.write("\n")
        post_p.write( " %8.5f" % (p[i]) )
        post_p.write("\n")

##########  End of Main Function  ##########        
############################################

rnd = Random()
sigma_q = 0.1
sigma_p = 0.0               
q_mean = -0.2	                      
p_mean = 0.0                 

#M = 1000                           
#q = []; p = []
#for i in range(M):
#    q.append(q_mean + sigma_q * rnd.normal() )
#    p.append(p_mean + sigma_p * rnd.normal() )
#    q.append(-1.35 + 0.005*i)
#    p.append(0.0)

q =  [-0.18453, -0.27176, -0.09767, -0.20181, -0.39299, -0.14552, -0.34007, -0.23812, -0.28699, 0.00468, -0.29117, -0.14674, -0.1673, -0.02548, -0.28093, -0.39157, -0.01229, -0.05718, -0.08138, -0.21451, -0.20842, -0.09283, -0.38671, -0.18711, -0.25665, -0.04732, -0.07632, -0.17876, -0.33484, -0.12246, -0.16006, -0.29642, -0.18408, -0.20119, -0.13925, -0.33054, -0.3183, -0.19843, -0.02648, -0.01383, -0.14895, -0.22899, -0.13109, -0.19454, -0.37925, -0.22662, -0.15791, -0.04113, -0.19493, -0.23624, -0.08741, -0.14928, -0.17918, -0.06081, -0.19778, -0.15588, -0.17809, -0.19984, -0.11012, -0.12469, 0.01132, -0.37729, -0.22897, -0.30538, -0.11951, -0.22946, -0.09377, -0.19522, -0.18005, -0.2536, 0.01374, -0.24069, -0.20997, -0.24385, -0.1052, -0.03915, -0.1484, -0.27224, -0.09733, -0.24751, -0.27639, -0.24795, -0.18712, -0.16843, -0.22368, -0.27533, -0.32934, -0.25385, -0.09588, -0.38511, -0.24775, 0.0324, -0.26857, -0.05366, -0.27866, -0.29449, -0.14433, 0.01426, -0.16376, -0.12976, -0.31253, -0.35133, -0.14692, -0.16377, -0.0459, -0.24752, -0.01625, 0.00136, -0.08057, -0.25746, -0.10331, -0.33956, -0.15575, -0.09775, -0.10801, -0.29143, -0.21169, -0.10895, -0.18949, -0.01341, -0.0842, -0.13775, -0.0805, -0.28336, -0.22546, -0.26765, 0.02211, 0.01285, -0.02677, -0.05918, -0.12651, -0.18813, -0.22532, -0.28469, -0.27949, -0.37248, -0.12167, -0.13301, -0.19049, -0.19369, -0.17276, -0.27361, -0.03607, -0.19393, -0.24393, -0.28116, -0.27569, -0.22222, -0.18494, -0.27616, -0.28723, -0.14514, -0.16355, -0.01938, -0.10415, -0.3204, -0.13695, -0.39696, -0.35615, -0.1106, -0.13398, -0.36393, -0.27267, -0.17289, -0.10515, -0.19349, -0.10041, -0.27616, -0.13383, -0.37622, -0.13306, -0.36372, -0.28087, -0.36581, -0.03764, -0.14551, -0.13669, -0.21585, -0.239, -0.07287, -0.10534, -0.29425, -0.31013, -0.15949, -0.18038, -0.19212, -0.18649, -0.25875, -0.24063, -0.12227, -0.29897, -0.24133, -0.36386, -0.15615, -0.24347, -0.33891, -0.26006, -0.34833, -0.24711, -0.22155, -0.37989, -0.05184, -0.34688, -0.26303, -0.23864, -0.26892, -0.09617, -0.25085, -0.15615, -0.16425, -0.26063, -0.20123, -0.12991, -0.14706, -0.24822, -0.42936, -0.17002, -0.07027, -0.17773, -0.1392, -0.22454, -0.03695, -0.20363, -0.00742, -0.14423, -0.15582, 0.02527, -0.14288, -0.36875, -0.17005, -0.37008, -0.24362, -0.12539, -0.21023, -0.09397, -0.12444, -0.13783, -0.19607, -0.23114, -0.32097, -0.13688, -0.14775, -0.05112, -0.13474, -0.20069, -0.24974, -0.18318, -0.36007, -0.19375, -0.1706, -0.0972, -0.1626, -0.34777, -0.20423, -0.26353, -0.21071, -0.1819, -0.37174, -0.02508, -0.21683, -0.24926, -0.03747, -0.28022, -0.07927, -0.32903, -0.11012, -0.18543, -0.34196, -0.16605, -0.2864, -0.27439, -0.4066, -0.21134, -0.19512, -0.27089, -0.10457, -0.20747, -0.06598, -0.18855, -0.23053, -0.1079, -0.29693, -0.22526, -0.16341, -0.23458, -0.1833, -0.15925, -0.25857, -0.2069, -0.04249, -0.34725, -0.25564, -0.30467, -0.25014, -0.1169, -0.14405, -0.18455, -0.24436, -0.08586, -0.26632, -0.30097, -0.35387, -0.19588, -0.23753, -0.27605, -0.24331, -0.34886, -0.09293, -0.0454, -0.17538, -0.2148, -0.18737, -0.30392, -0.37271, -0.23819, -0.28199, -0.10057, -0.22604, -0.2456, -0.11972, -0.0969, -0.14093, -0.31158, -0.08876, -0.25464, -0.23466, -0.29193, -0.09695, -0.27806, -0.21694, -0.05689, -0.26239, -0.22437, -0.09202, -0.16327, -0.18092, -0.12589, -0.13931, -0.24474, -0.28705, -0.27784, -0.29495, -0.2125, -0.21679, -0.23531, -0.14022, -0.30853, -0.23944, -0.24226, -0.28264, -0.0217, -0.15382, -0.33185, -0.26534, -0.2934, -0.18915, -0.3461, -0.13539, -0.21631, -0.21945, -0.24532, -0.22091, -0.28098, 0.00859, -0.15567, -0.06592, -0.2459, -0.09308, -0.1759, -0.1277, -0.21483, -0.29595, -0.05126, -0.19801, -0.1015, -0.10336, -0.26091, -0.26626, -0.11824, -0.05968, -0.40093, -0.20932, -0.23216, -0.17383, -0.18219, -0.01334, -0.13451, -0.26189, -0.06823, -0.44875, -0.31005, -0.23165, -0.2629, -0.2266, -0.20967, -0.33711, -0.34647, -0.19268, -0.24991, -0.26537, -0.33243, -0.31942, -0.14145, -0.06123, -0.2514, -0.09546, -0.03618, -0.20724, -0.14539, -0.23606, -0.06113, -0.10447, -0.2334, -0.14272, -0.18946, -0.17433, -0.11989, -0.26204, -0.07714, -0.12458, -0.28294, -0.0804, -0.21006, -0.20153, -0.26167, -0.13911, -0.26112, -0.12268, -0.25251, -0.20357, -0.13234, -0.112, -0.25899, -0.20939, -0.22834, -0.32716, -0.16895, -0.2321, -0.13433, -0.33444, -0.12134, -0.10989, -0.37321, -0.40063, -0.12926, -0.37645, -0.06648, -0.26751, -0.44619, -0.29473, -0.30642, -0.23606, -0.23135, -0.26232, -0.16171, -0.15314, -0.19231, -0.22855, -0.2305, -0.41389, -0.25637, -0.06871, -0.41909, -0.2184, -0.23585, -0.25646, -0.07829, -0.2911, -0.1403, -0.09292, -0.15069, -0.28733, -0.40178, -0.41283, -0.38425, -0.04764, -0.23542, -0.27985, 0.01621, 0.03653, -0.24979, -0.20378, -0.39135, -0.28339, -0.23481, -0.23618, -0.15245, -0.25053, -0.11842, -0.10356, -0.17917, -0.18381, -0.35121, -0.23915, -0.04527, -0.26219, -0.34731, -0.07504, -0.18551, -0.15446, -0.16459, -0.23397, -0.41503, -0.16233, -0.20016, -0.19271, -0.29592, -0.31316, -0.27908, -0.28962, -0.18177, -0.20828, -0.15913, -0.24216, -0.03945, -0.16048, -0.14788, -0.38215, -0.21457, -0.15888, -0.20846, -0.29309, 0.08447, -0.22863, -0.01254, -0.15021, -0.24731, -0.17355, -0.33844, -0.04716, -0.18547, 0.01497, -0.19091, -0.17442, -0.12082, -0.08164, -0.05508, -0.06227, -0.26288, -0.3279, -0.06658, -0.32281, -0.35726, -0.33909, -0.34581, -0.31005, -0.18341, -0.16155, -0.28676, -0.30161, -0.22466, -0.23982, -0.17326, -0.13177, -0.37976, -0.28224, -0.23402, -0.10604, -0.19239, -0.38329, -0.26154, -0.18129, -0.18077, -0.29968, -0.28403, -0.08045, -0.2002, -0.25419, -0.22159, -0.13212, -0.17165, -0.1609, -0.32192, -0.331, -0.25551, -0.11217, -0.26156, -0.19594, -0.17079, -0.24534, -0.06975, -0.17337, -0.17161, -0.24311, -0.03761, -0.26029, -0.3994, -0.27483, -0.2363, -0.22838, -0.11808, -0.28048, -0.20767, -0.253, -0.08618, -0.18985, -0.18461, -0.10401, -0.24241, 0.00985, -0.29056, -0.33222, -0.1921, -0.23994, -0.33804, -0.14086, -0.15158, -0.09543, -0.15681, -0.09667, -0.31385, -0.22022, -0.17624, -0.17395, -0.0938, -0.20585, -0.04455, -0.18453, -0.12454, -0.164, -0.35203, -0.16799, -0.17056, -0.14951, -0.07032, 0.01692, -0.20662, -0.139, -0.20585, -0.09696, 0.01257, -0.17856, -0.22378, -0.28907, -0.1813, -0.1552, -0.19551, -0.20278, -0.27246, -0.30839, -0.26843, -0.30566, -0.15007, -0.03922, -0.10968, -0.24167, -0.23642, -0.30859, -0.22103, -0.08988, -0.05742, -0.24682, -0.12714, -0.33029, -0.18119, -0.22648, -0.23694, -0.22031, -0.07219, -0.11182, -0.24119, -0.20337, -0.09851, -0.39392, -0.40708, -0.27605, -0.19133, -0.26516, -0.27003, -0.33111, -0.12452, -0.2224, -0.30412, -0.13237, -0.22476, -0.2797, -0.07882, -0.3418, -0.17556, -0.1625, -0.09734, -0.22368, -0.07121, 0.01663, -0.32997, -0.17467, -0.24702, -0.17869, -0.2449, -0.11426, -0.22216, -0.01272, -0.25191, -0.0671, -0.17352, -0.24571, -0.01589, -0.39609, -0.28963, -0.03655, -0.31173, -0.30335, -0.06037, -0.13964, -0.2181, -0.21818, -0.22678, -0.43345, -0.1519, -0.10976, -0.27989, -0.36617, -0.17293, -0.20544, -0.39469, -0.26352, -0.1347, -0.11777, -0.41231, -0.34826, -0.12063, -0.26781, -0.13306, -0.20221, -0.20984, -0.00455, -0.19549, -0.13362, -0.3962, -0.12501, -0.3359, -0.16727, -0.37211, -0.05469, -0.17151, -0.31654, -0.18137, -0.13335, -0.22945, -0.23432, -0.2507, -0.22617, -0.11569, -0.16843, -0.28591, -0.12185, -0.29151, -0.40593, -0.30739, -0.31028, -0.14845, -0.06193, -0.15567, -0.20661, -0.16903, -0.31671, -0.09806, -0.26882, -0.15472, -0.20733, -0.2847, -0.39849, -0.27295, -0.16914, -0.17259, -0.20746, -0.29202, 0.03028, -0.38484, -0.21282, -0.2454, -0.27216, -0.31859, 0.06705, -0.22306, -0.28206, -0.15723, -0.23268, -0.24182, -0.16151, -0.23632, -0.15863, -0.119, -0.185, 0.00391, -0.34749, -0.14798, -0.14131, -0.22253, -0.08208, -0.15577, -0.25101, -0.24222, -0.21458, -0.33798, -0.21058, -0.24754, -0.52511, -0.26443, -0.03602, 0.00625, -0.31389, -0.05266, -0.20267, -0.10655, -0.10437, -0.12775, -0.23841, -0.22411, -0.14995, -0.2147, -0.29628, -0.3196, -0.20141, -0.21003, -0.25529, -0.21204, -0.13701, -0.16748, -0.21172, -0.23633, -0.28915, -0.23922, -0.18383, -0.11997, -0.01384, -0.05829, -0.39373, -0.26966, 0.00595, -0.20552, -0.32952, -0.13494, -0.16915, 0.03917, -0.08644, -0.08341, -0.16385, -0.21967, -0.39882, -0.1335, -0.21544, -0.1853, -0.19125, -0.28557, -0.22401, -0.34015, -0.25631, -0.43011, -0.0523, -0.44493, -0.21767, -0.0273, -0.17847, -0.11817, -0.18326, -0.28712, -0.12226, -0.15046, -0.19608, -0.1599, -0.31506, 0.04411, -0.01524, -0.13501, -0.20614, -0.34586, -0.14774, -0.20707, -0.21929, -0.15391, -0.21752, -0.06828, -0.42704, -0.0955, -0.38335, -0.19262, -0.05885, -0.24404, -0.38056, -0.09788, -0.15276, -0.27607, -0.41071, -0.14226, -0.25394, -0.20083, -0.29679, -0.05201, -0.29032, -0.22879, -0.06414, -0.0352, -0.20251, -0.25038, -0.26943, -0.24085, -0.04956, -0.19828, -0.03, 0.01668, -0.17478, -0.34188, -0.18755, -0.10278, -0.12484, -0.34083, -0.22915, -0.32544, -0.12652, -0.28274, -0.15686, -0.21345, -0.11984, -0.22808, -0.146, -0.17726, -0.39845, -0.23514, -0.08484, -0.1324, -0.12784, -0.16545, -0.39211, -0.1009, -0.08535, -0.14868, -0.16745, -0.27333, -0.22774, -0.0421, -0.12164, -0.21191, -0.13507, -0.18955, -0.08228, -0.23257, -0.13948, -0.1746, -0.20771, -0.14663, -0.44111, -0.16971, -0.08937, -0.12893, -0.10619, -0.12902, -0.27464, -0.23356, -0.29661, -0.15484, -0.31087, -0.11882, -0.25295, -0.5606, -0.29346, -0.11397, -0.1397, -0.18055, -0.21407, -0.09022, -0.19472, -0.14299, -0.17273, -0.20038, -0.22744, -0.18046, -0.13825, -0.19861, -0.18016, -0.04925, -0.23774, -0.28286, -0.25381, -0.27819, -0.3523, -0.25169, -0.16525, -0.26595, -0.18962, -0.1825, -0.02164, -0.1578, -0.27605, -0.12276, -0.17588, -0.01071, -0.28836, -0.17415, -0.10834, 0.00052, -0.28108, -0.15315, -0.15549, -0.25247, -0.42705, -0.26649, -0.42155, -0.10107, -0.12583]

p =  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

q_grid = []; p_grid = []
for i in range(-150,200):
    q_grid.append(0.01*i)

for i in range(-100,500):
    p_grid.append(0.1*i)                                

ent_type = 1      # 1 - ETHD, 2 - RPMD, 3 - QHD2
model = 1         # 1 - cubic, 2 - double well 
vel_rescale = 0   # 0 - no, 1 - yes
approx = 2        # 1 - H = H0,  2 - H = H0 + H1
dt = 0.1
m = 2000.0
Nsnap = 50000
Nstep = 1
#                          
main(q, p, q_grid, p_grid, Nstep, Nsnap, m, dt, approx, model, ent_type, vel_rescale)  
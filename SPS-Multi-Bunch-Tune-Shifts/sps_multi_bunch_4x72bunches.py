from __future__ import division

import sys, os
BIN = "/hpcscratch/user/imasesso/PyHEADTAIL_py3/PyHEADTAIL/" # Path that contains PyHEADTAIL directory
sys.path.append(BIN)

import mpi4py.rc
mpi4py.rc.threaded = False

import time
import numpy as np
from mpi4py import MPI
from scipy.constants import c, e, m_p

from PyHEADTAIL.monitors.monitors import BunchMonitor, SliceMonitor, ParticleMonitor
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from PyHEADTAIL.impedances.wakes import WakeField, CircularResonator, Resonator, WakeTable
from PyHEADTAIL.feedback.feedback import IdealBunchFeedback
from PyHEADTAIL.machines.synchrotron import Synchrotron
from SPSoctupolesNewConfiguration import SPSOctupolesNew
import PyHEADTAIL.aperture.aperture as aperture

# The parameters in this simulation correspond to SPS injection energy (26Gev) for Q-20 optics.
# Simulation example for LHC-type beams in the SPS

def run(xi):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    print('I am rank ' + str(rank) + ' out of ' + str(size))
    print('  ')
    np.random.seed(0+rank)

    # BEAM AND MACHINE PARAMETERS
    # ============================
    intensity = xi # *bunch* intensity
    n_turns = 500
    n_macroparticles = 8e5 # number of macroparticles *per bunch*

    charge = e
    mass = m_p
    alpha = 0.0030864197530864196

    gamma = 27.728550064942279
    p0 = np.sqrt(gamma**2 - 1.) * mass * c

    accQ_x = 20.10
    accQ_y = 20.28
    #Q_s = 0.017046597789257857
    circumference = 6911.5038378975451

    s = None
    alpha_x = None
    alpha_y = None
    beta_x = 54.644808743169399
    beta_y = 54.509415262636274
    D_x = 0
    D_y = 0
    optics_mode = 'smooth'
    name = None
    n_segments = 1

    # Carlo's implementation of octupoles
    KLOF=0
    KLOD=0
    Octupoles=SPSOctupolesNew()
    app_x =2 * p0 * Octupoles.get_anharmonicities(KLOF,KLOD)[0]
    app_xy=2 * p0 * Octupoles.get_anharmonicities(KLOF,KLOD)[1]
    app_y =2 * p0 * Octupoles.get_anharmonicities(KLOF,KLOD)[2]
    Qp_x = np.array([0.4*accQ_x, Octupoles.get_q2(KLOF,KLOD)[0]+272,-1869000]) 
    Qp_y = np.array([0.45*accQ_y, Octupoles.get_q2(KLOF,KLOD)[1]+662,1449600])

    # Aperture
    apt_xy = aperture.EllipticalApertureXY(x_aper=0.07, y_aper=0.03)

    # RF parameters 
    longitudinal_mode = 'non-linear'
    h_RF = [ 4620, 4*4620 ]
    V_RF = [ 4.5e6, 0.45e6 ]
    dphi_RF = [ 0., np.pi ]
    wrap_z = False

    machine = Synchrotron(
            optics_mode=optics_mode, circumference=circumference,
            n_segments=n_segments, s=s, name=name,
            alpha_x=alpha_x, beta_x=beta_x, D_x=D_x,
            alpha_y=alpha_y, beta_y=beta_y, D_y=D_y,
            accQ_x=accQ_x, accQ_y=accQ_y, Qp_x=Qp_x, Qp_y=Qp_y,
            app_x=app_x, app_y=app_y, app_xy=app_xy,
            alpha_mom_compaction=alpha, longitudinal_mode=longitudinal_mode,
            h_RF=np.atleast_1d(h_RF), V_RF=V_RF, dphi_RF=dphi_RF, p0=p0, p_increment=0.,
            charge=charge, mass=mass, wrap_z=wrap_z)#, Q_s=Q_s)

    h = np.min(machine.longitudinal_map.harmonics) * 1

    # FILLING SCHEME
    # ==============
    h_bunch = 924 # does not necessarily correspond to h (h_rf), but has to correspond to smallest bunch
                  # separation. Here, 4620 / 924 = 5, hence 5 rf buckets. See also bunch_spacing variable.
                  # Important for performance with FFT algorithms, since whole ring is sliced.
    n_batches = 4 # Number of batches
    batch_length = 72 # Number of bunches per batch
    bunch_spacing = 1 # Gap between bunches in units of harmonic cycles, i.e. h_bunch
    batch_gap = 10 # Gap between batches in units of bunch spacing
    offset = 0 # Location of the first bunch in units of bunch spacing

    filling_scheme = []
    for j in range(n_batches):
        for i in range(batch_length):
            filling_scheme.append(bunch_spacing*(offset + i+j*(batch_length+batch_gap)))

    # BEAM
    # ====
    epsn_x = 1.8e-6
    epsn_y = 1.8e-6
    sigma_z = 0.23
    allbunches = machine.generate_6D_Gaussian_bunch(
        n_macroparticles, intensity, epsn_x, epsn_y, sigma_z=sigma_z,
        filling_scheme=filling_scheme, matched=True)	#was False in linear mode
   

    # CREATE BEAM SLICERS
    # ===================
    # IMPORTANT TO USE FIXED LONGITUDINAL CUTS!
    # Pass correct h_bunch, also considered above when defining filling scheme
    slicer_for_wakefields = UniformBinSlicer(300, z_cuts=(-3.*sigma_z, 3.*sigma_z),
                                             circumference=circumference, h_bunch=h_bunch)

    # WAKEFIELD
    # =========
    mpi_settings = 'linear_mpi_full_ring_fft'

    # n_turns_wake: memory for wake field. 
    n_turns_wake = 43 

    wakefile1 = ('/hpcscratch/user/imasesso/simulations/wakes/SPS_wake_model_with_EF_BPH_BPV_kicker_wall_steps_ZS_PostLS2_Q20.txt')
    ww1 = WakeTable(wakefile1, ['time', 'dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y'], n_turns_wake=n_turns_wake)
    wake_field = WakeField(slicer_for_wakefields, ww1, mpi=mpi_settings,
                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
    machine.one_turn_map.append(wake_field)


    # TRANSVERSE FEEDBACK SYSTEM
    # ==========================
    dampingrate = 32.
    damper = IdealBunchFeedback(gain=2./dampingrate, multi_bunch=True)
    machine.one_turn_map.append(damper)


    # MONITORS
    # ========
    # There is one bunchmonitor file for the whole beam.
    # f = h5py.File('bunchmonitor.h5')
    # f['Bunches'].keys() will give the keys for all the bunches that are stored in file.
    # To plot mean_x of "3rd bunch" for example: plt.plot(f['Bunches']['3']['mean_x'])

    write_every=100
    bunchmonitor = BunchMonitor(
         './bunchmonitor_coherent_tunes_{:d}wake_{:d}batches_{:d}bunches_xi{:.0f}'.format(int(n_turns_wake),int(n_batches),int(batch_length),xi), n_turns, simulation_parameters_dict={}, write_buffer_every=write_every, buffer_size=n_turns,mpi=True, filling_scheme=filling_scheme)

 
    # TRACKING LOOP
    # =============
    print('\n--> Begin tracking...\n')

    for i in range(n_turns):
        t0 = time.perf_counter()

        if i == 172: # To kick all bunches at a certain turn number in the horizontal plane 
            allbunches.x += 1e-4
        if i == 301: # To kick all bunches at a certain turn number in the vertical plane 
            allbunches.y += 1e-4

        machine.track(allbunches)
        apt_xy.track(allbunches)

        # Dump bunch and slice monitor data
        bunchmonitor.dump(allbunches)

        if i%20 != 0:
            continue

        if rank == 0:
            print('Turn: {:4d}, \tTime: {:3s}'.format(i, str(time.perf_counter() - t0)))
    print('\n*** Successfully completed!')


if __name__=="__main__":
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    xi_vect = np.array([1.8e11, 2.0e11, 2.2e11, 2.2e11])

    for xi in xi_vect:
        run(xi=xi)

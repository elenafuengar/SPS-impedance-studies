

import sys, os

import sys, time
import numpy as np
from scipy.constants import c, e, m_p

from PyHEADTAIL.monitors.monitors import BunchMonitor, SliceMonitor, ParticleMonitor
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from PyHEADTAIL.impedances.wakes import WakeField, CircularResonator, Resonator, WakeTable
from PyHEADTAIL.feedback.transverse_damper import TransverseDamper
from PyHEADTAIL.machines.synchrotron import Synchrotron
#from SPSoctupolesNewConfiguration import SPSOctupolesNew
import PyHEADTAIL.aperture.aperture as aperture
#apt_xy = aperture.EllipticalApertureXY(x_aper=0.07, y_aper=0.03)
#self.machine.one_turn_map.append(apt_xy)

from tqdm import tqdm

def run(xi,target_tune,n_slices):
    # BEAM AND MACHINE PARAMETERS
    # ============================
    intensity = xi # *bunch* intensity
    n_turns = 500
    n_macroparticles = 300000 # number of macroparticles *per bunch*

    charge = e
    mass = m_p
    alpha = 0.0030864197530864196

    gamma = 27.728550064942279
    p0 = np.sqrt(gamma**2 - 1.) * mass * c

    accQ_x = 20.15
    accQ_y = target_tune
    fracQ_y = int(10000*(np.modf(accQ_y)[0]))
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

    # Chromaticity and amplitude detuning from octupoles
#    Qp_x = xi*accQ_x
#    Qp_y = xi*accQ_y

# MS implementation
#    axx=-500
#    ayy=-1000
#    app_x = 2. * p0 * axx
#    app_y = 2. * p0 * ayy
#    app_xy = 0.

# My implementation of octupoles
    KLOF=0
    KLOD=0
    #Octupoles=SPSOctupolesNew()
    app_x = 0 #2 * p0 * Octupoles.get_anharmonicities(KLOF,KLOD)[0]
    app_xy= 0 #2 * p0 * Octupoles.get_anharmonicities(KLOF,KLOD)[1]
    app_y = 0 #2 * p0 * Octupoles.get_anharmonicities(KLOF,KLOD)[2]
    Qp_x = 0 #np.array([0.4*accQ_x, Octupoles.get_q2(KLOF,KLOD)[0]+272,-1869000])
    Qp_y = 0 #np.array([0.4*accQ_y, Octupoles.get_q2(KLOF,KLOD)[1]+662,1449600])

#    Qp_x = np.array([xi*accQ_x+Octupoles.get_q1_feeddown(KLOF,KLOD)[0],Octupoles.get_q2(KLOF,KLOD)[0]]) Include also effect on first order chroma
#    Qp_y = np.array([xi*accQ_x+Octupoles.get_q1_feeddown(KLOF,KLOD)[1],Octupoles.get_q2(KLOF,KLOD)[1]])

    apt_xy = aperture.EllipticalApertureXY(x_aper=0.07, y_aper=0.03)
    # self.machine.one_turn_map.append(apt_xy)


    longitudinal_mode = 'non-linear'

    # RF parameters (not really relevant here, since linear tracking)
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

    # BEAM
    # ====
    epsn_x = 1.5e-6
    epsn_y = 1.5e-6
    sigma_z = 0.23
    bunch = machine.generate_6D_Gaussian_bunch(
        n_macroparticles, intensity, epsn_x, epsn_y, sigma_z=sigma_z)	#was False in linear mode
    bunch.y +=1e-3
    bunch.x +=1e-3
    # CREATE BEAM SLICERS
    # ===================
    # IMPORTANT TO USE FIXED LONGITUDINAL CUTS!
    # Pass correct h_bunch, also considered above when defining filling scheme
    h_bunch = 924
    slicer_for_wakefields = UniformBinSlicer(n_slices, z_cuts=(-3.*sigma_z, 3.*sigma_z),
                                             circumference=circumference, h_bunch=h_bunch)	#it was 40 and 20 and h_bunch=h_bunch
    slicer_for_monitor = UniformBinSlicer(80, z_cuts=(-4.*sigma_z, 4.*sigma_z),
                                          circumference=circumference, h_bunch=h_bunch)


    # WAKEFIELD
    # =========
    # Depending on number of bunches in beam: use 'memory_optimized', 'linear_mpi_full_ring_fft',
    # (or 'circular_mpi_full_ring_fft') for best performance (see slides).
    # Circular_fft not robust when using quad wakes, or when other elements in ring.
    # General advice: use either memory_optimized or linear_mpi_full_ring_fft.

    # n_turns_wake: memory for wake field. See slides for convention (has changed compared to HEADTAIL).
    n_turns_wake = 43 

 # BBR -- for kickers
#    bbr_h_dip = Resonator(R_shunt=1.6e6, frequency=1.3e9, Q=1., Yokoya_X1=1., Yokoya_Y1=0., Yokoya_X2=0., Yokoya_Y2=0., switch_Z=False, n_turns_wake=n_turns_wake)
#    bbr_h_quad = Resonator(R_shunt=1.6e6*1.0, frequency=0.6e9, Q=1., Yokoya_X1=0., Yokoya_Y1=0., Yokoya_X2=-1., Yokoya_Y2=0., switch_Z=False, n_turns_wake=n_turns_wake)
#    bbr_v_dip = Resonator(R_shunt=2.3e6, frequency=0.9e9, Q=1., Yokoya_X1=0., Yokoya_Y1=1., Yokoya_X2=0., Yokoya_Y2=0., switch_Z=False, n_turns_wake=n_turns_wake)
#    bbr_v_quad = Resonator(R_shunt=1.7e6*1.0, frequency=0.6e9, Q=1., Yokoya_X1=0., Yokoya_Y1=0., Yokoya_X2=0., Yokoya_Y2=1., switch_Z=False, n_turns_wake=n_turns_wake)

    # NBR -- for kickers
#    nbr_h_dip = Resonator(R_shunt=7.3e6, frequency=74.5e6, Q=4., Yokoya_X1=1., Yokoya_Y1=0., Yokoya_X2=0., Yokoya_Y2=0., switch_Z=False, n_turns_wake=n_turns_wake)
#    nbr_h_quad = Resonator(R_shunt=14e6*1.0, frequency=43.0e6, Q=4., Yokoya_X1=0., Yokoya_Y1=0., Yokoya_X2=-1., Yokoya_Y2=0., switch_Z=False, n_turns_wake=n_turns_wake)
#    nbr_v_dip = Resonator(R_shunt=5.7e6, frequency=80.0e6, Q=4., Yokoya_X1=0., Yokoya_Y1=1., Yokoya_X2=0., Yokoya_Y2=0., switch_Z=False, n_turns_wake=n_turns_wake)
#    nbr_v_quad = Resonator(R_shunt=11.0e6*1.0, frequency=46.0e6, Q=5., Yokoya_X1=0., Yokoya_Y1=0., Yokoya_X2=0., Yokoya_Y2=1., switch_Z=False, n_turns_wake=n_turns_wake)


#    nbr_RF200 = CircularResonator(R_shunt=155e6, frequency=939.45e6, Q=8500., n_turns_wake=n_turns_wake)
    '''
    wakefile1 = ('/hpcscratch/user/imasesso/simulations/wakes/kickerSPSwake_2020_MKP_MB.wake')
    wakefile2 = ('/hpcscratch/user/imasesso/simulations/wakes/WallWake_CZ_Q20_SPS.wake')
#    wakefile3 = ('/hpcscratch/user/imasesso/simulations/wakes/MB_step_wake_analytical.wake')
    ww1 = WakeTable(wakefile1, ['time', 'dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y'], n_turns_wake=n_turns_wake)
    ww2 = WakeTable(wakefile2, ['time', 'dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y'], n_turns_wake=n_turns_wake)
#    ww3 = WakeTable(wakefile3, ['time', 'dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y'], n_turns_wake=n_turns_wake)

    #wake_source_list = [bbr_kicker, ww]
    #wake_field = WakeField(slicer_for_wakefields, wake_source_list, mpi=mpi_settings,
    #                       Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
    wake_field_kicker = WakeField(slicer_for_wakefields, ww1, mpi=mpi_settings,
                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)

    wake_field_wall = WakeField(slicer_for_wakefields, ww2, mpi=mpi_settings,
                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)

#    wake_field_step = WakeField(slicer_for_wakefields, ww3, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)

#    wake_field_bbrhdip = WakeField(slicer_for_wakefields, bbr_h_dip, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    wake_field_bbrhquad = WakeField(slicer_for_wakefields, bbr_h_quad, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    wake_field_bbrvdip = WakeField(slicer_for_wakefields, bbr_v_dip, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    wake_field_bbrvquad = WakeField(slicer_for_wakefields, bbr_v_quad, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    wake_field_nbrhdip = WakeField(slicer_for_wakefields, nbr_h_dip, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    wake_field_nbrhquad = WakeField(slicer_for_wakefields, nbr_h_quad, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    wake_field_nbrvdip = WakeField(slicer_for_wakefields, nbr_v_dip, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    wake_field_nbrvquad = WakeField(slicer_for_wakefields, nbr_v_quad, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)

    machine.one_turn_map.append(wake_field_kicker)
    machine.one_turn_map.append(wake_field_wall)
    '''
#    machine.one_turn_map.append(wake_field_step)
#    machine.one_turn_map.append(wake_field_bbrhdip)
#    machine.one_turn_map.append(wake_field_bbrhquad)
#    machine.one_turn_map.append(wake_field_bbrvdip)
#    machine.one_turn_map.append(wake_field_bbrvquad)
#    machine.one_turn_map.append(wake_field_nbrhdip)
#    machine.one_turn_map.append(wake_field_nbrhquad)
#    machine.one_turn_map.append(wake_field_nbrvdip)
#    machine.one_turn_map.append(wake_field_nbrvquad)


#    wake_field_RF200 = WakeField(slicer_for_wakefields, nbr_RF200, mpi=mpi_settings,
#                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)
#    machine.one_turn_map.append(wake_field_full)
#    machine.one_turn_map.append(wake_field_RF200)
#    machine.one_turn_map.append(apt_xy)

    # TRANSVERSE FEEDBACK SYSTEM
    # ==========================
    dampingrate = 30.
    #damper = IdealBunchFeedback(gain=2./dampingrate, multi_bunch=True)
    damper = TransverseDamper(dampingrate_x=dampingrate, dampingrate_y=dampingrate)
    machine.one_turn_map.append(damper)


    # MONITORS
    # ========
    # There is one bunchmonitor file for the whole beam.
    # f = h5py.File('bunchmonitor.h5')
    # f['Bunches'].keys() will give the keys for all the bunches that are stored in file.
    # To plot mean_x of "3rd bunch" for example: plt.plot(f['Bunches']['3']['mean_x'])
    write_every=100
    n_batches=1
    batch_length=1
    bunchmonitor = BunchMonitor(
         './bunchmonitor_{:d}wake_{:d}batch_{:d}bunches_{:d}setTune_{:d}Slices_xi{:.0f}'.format(int(n_turns_wake),int(n_batches),int(batch_length), fracQ_y,n_slices,xi), n_turns,
         simulation_parameters_dict={}, write_buffer_every=n_turns, buffer_size=n_turns)

    # Slicemonitor cannot yet handle multi-bunch beams. Need to do manully with split_to_views()
    # as shown below. 1 file per bunch is created, but can also select smaller number of bunches.
    slicemonitor_dict = {}
    n_trns_slice_max = 50
    bid = 0 #bunch.bucket_id[0]
    slicemonitor_dict[str(bid)] = SliceMonitor(
            './slicemonitor_{:d}wake_{:d}batch_{:d}bunches_{:d}setTune_{:d}Slices_b{:d}_xi{:.1f}'.format(int(n_turns_wake),int(n_batches),int(batch_length),fracQ_y,n_slices,int(bid),xi), n_trns_slice_max, slicer_for_monitor)

#    particlemonitor_dict = {}    # Only save 1every tenth particle ... stride defines the stepsize (i.e. how many
   				 # particles are skipped when saving to monitor)
#    stride=1
#    bunch_list = allbunches.split_to_views()
#    for b in bunch_list[10:11]:
#        bid = b.bucket_id[0]
#        particlemonitor_dict[str(bid)] = ParticleMonitor('./particlemonitor_b{:d}_xi{:.1f}'.format(int(bid), xi), stride=stride)

    # TRACKING LOOP
    # =============
    print('\n--> Begin tracking...\n')
    n_trns_slice = 0
    activated_slicemonitor = False

    for i in tqdm(list(range(n_turns))):
        t0 = time.time()

        machine.track(bunch)
        # apt_xy.track(allbunches)

        # Dump bunch and slice monitor data
        bunchmonitor.dump(bunch)

        # This is only for slice monitor 
        apt_xy.track(bunch)

        if activated_slicemonitor == False:
            
            if ((bunch.mean_x() > 1.) or (n_trns_slice > 0) or i > (n_turns-n_trns_slice_max)):
                activated_slicemonitor = True

            if activated_slicemonitor and (n_trns_slice < n_trns_slice_max):
                slicemonitor_dict[str(0)].dump(bunch)
                n_trns_slice += 1

#        bunch_list = allbunches.split_to_views()
#        for jj, bb in enumerate(bunch_list[10:11]):
#            particlemonitor_dict[str(bb.bucket_id[0])].dump(bb)

        #if i%20 != 0:
        #    continue

        #print(('Turn: {:4d}, \tTime: {:3s}'.format(i, str(time.time() - t0))))
    print('\n*** Successfully completed!')


if __name__=="__main__":
    xi_vect    = np.array([0.2e11, 0.3e11, 0.4e11,0.5e11, 0.6e11, 0.7e11, 0.8e11, 0.9e11, 1.0e11, 1.1e11, 1.2e11, 1.3e11, 1.4e11, 1.5e11, 1.6e11, 1.7e11, 1.8e11, 1.9e11, 2.0e11, 2.1e11, 2.2e11, 2.3e11, 2.4e11, 2.5e11, 2.6e11])
    #xi_vect = [1.4e11]
    #tunes_vect = np.array([20.202, 20.212, 20.223, 20.232, 20.241, 20.250, 20.257, 20.264, 20.271, 20.277, 20.283, 20.288, 20.293])
    #xi_vect = np.array([2.2e11, 2.4e11, 2.6e11])
    target_tune = 20.40
    n_slices_vect = [50]
    for n_slices in n_slices_vect:
        for xi in xi_vect:
            run(xi=xi, target_tune=target_tune,n_slices=n_slices)

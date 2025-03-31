import numpy as np
import h5py
import os
from tqdm import tqdm
import glob, shutil

from scipy.constants import c, e,  m_p
from SPSoctupolesNewConfiguration import SPSOctupolesNew

import xtrack as xt
import xobjects as xo
import xpart as xp
import xfields as xf




def run(qpv_val=0, n_slices=300, n_mp=1e6, intensity=3e10):
# Simulation settings
    n_turns = 6000
    n_turns_wake = 5
    n_macroparticles = int(n_mp)
    n_slices_wakes = 300
    average_values = True

    circumference = 6911.5038378975451

    averageRadius = circumference / (2 * np.pi)

    gamma=27.64
    beta = np.sqrt(1-1/gamma**2)
    f_rev = beta*c/circumference    #freq_sps = 43.2e3

    bl_val_ns=3.2
    sigma_z = bl_val_ns*1.e-9*c/4

    bucket_length_m = 12*sigma_z

    p0c = 26e9
    p0=p0c * e / c

    accQx=26.13
    accQy=26.18

    beta_x = averageRadius/accQx
    beta_y = averageRadius/accQy

    qph=0.1
    qpv=qpv_val



    

    #Adding wakes from file

    wake_table_filename = 'new_model_corr.wake'
    
    wake_file_columns = ['time', 'dipole_x', 'dipole_y',
                        'quadrupole_x', 'quadrupole_y']
    waketable=xf.Wakefield.table_from_headtail_file( wake_table_filename, wake_file_columns, flag_pyht_units=True)
    waketable['time']*= -1
    waketable['dipole_x']*= -1
    waketable['dipole_y']*= -1
    waketable['quadrupole_x']*= -1
    waketable['quadrupole_y']*= -1

    use_components=['dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y']
    wakefield = xf.Wakefield.from_table(
        waketable,use_components,
        #use_components=['dipole_x', 'dipole_y','quadrupole_x', 'quadrupole_y'],
        zeta_range=(-0.5*bucket_length_m, 0.5*bucket_length_m),
        num_slices=n_slices_wakes,
        num_turns=5.,
        circumference=circumference
        )


    #Optical parameters

    Octupoles=SPSOctupolesNew(optics='Q26')
    app_x= 2 * p0 * Octupoles.get_anharmonicities(0,0)[0]
    app_xy= 2 * p0 * Octupoles.get_anharmonicities(0,0)[1]
    app_y= 2 * p0 * Octupoles.get_anharmonicities(0,0)[2]


    Qpx = np.array([accQx, qph*accQx, Octupoles.get_q2(0,0)[0]+612.213*2,-90760.404*6])
    Qpy = np.array([accQy, qpv*accQy , Octupoles.get_q2(0,0)[1]-55.902*2,+55572.440*6])


    #RF parameters

    h_RF = np.array([ 4620, 4*4620 ]) 
    f_RF = h_RF * f_rev
    V_RF = [4.5e6, 0.3e6]
    dphi_RF = [180 , 0 ]

    #building one turn map
    one_turn_map = xt.LineSegmentMap(
        length=circumference, betx=beta_x, bety=beta_y,
        #qx=accQx, qy=accQy,
        #dqx=qph, dqy=qpv,
        longitudinal_mode='nonlinear',
        momentum_compaction_factor=0.0019236688211757462,
        voltage_rf=V_RF, frequency_rf = f_RF,
        lag_rf = dphi_RF,
        
        det_xx = app_x/p0, det_xy = app_xy/p0, det_yy=app_y/p0,
        dnqx=Qpx, dnqy=Qpy

    )
    '''
    one_turn_map2 = xt.LineSegmentMap(
        length=circumference, betx=42., bety=42.,
        qx=26.13, qy=26.18,
        longitudinal_mode='linear_fixed_qs',
        dqx=0.1, dqy=-1.35, 
        qs=2e-3, bets=731.27
    )
    '''
    # Generate line
    line = xt.Line(elements=[one_turn_map, wakefield],
                element_names=['one_turn_map', 'wf'])

    particle_ref = xt.Particles(p0c=p0c, mass0=xt.PROTON_MASS_EV)
    line.particle_ref = particle_ref
    line.build_tracker()

    # Generate particles
    particles = xp.generate_matched_gaussian_bunch(line=line,
                        num_particles=n_macroparticles , total_intensity_particles=intensity,
                        nemitt_x=1.5e-6, nemitt_y=1.5e-6, sigma_z=sigma_z, particle_ref=particle_ref)






    filename = f'sps_q26_qpv{qpv}_avr.h5'
    HT_format = True
    HT_format_dict={
        'Bunch': {
            'mean_x': [],
            'mean_y': [],
            'mean_xp': [],
            'mean_yp': [],
        },

    }
    for i in tqdm(range(n_turns), desc ="OTM tracking"):
        
        if i==1:
            particles.y += 1e-3
        
        line.track(particles, num_turns=1)
        if HT_format:
            HT_format_dict['Bunch']['mean_x']=np.append(HT_format_dict['Bunch']['mean_x'],np.mean(particles.x))
            HT_format_dict['Bunch']['mean_y']=np.append(HT_format_dict['Bunch']['mean_y'],np.mean(particles.y))
            HT_format_dict['Bunch']['mean_xp']=np.append(HT_format_dict['Bunch']['mean_xp'],np.mean(particles.px))           
            HT_format_dict['Bunch']['mean_yp']=np.append(HT_format_dict['Bunch']['mean_yp'],np.mean(particles.py))

        else:
            monitor_to_h5(particles,filename,collapse_into_avr=True)
    if HT_format:
        with h5py.File(f'q26_mean_vals_qpv{qpv}_N{intensity}__nonlinear.h5', 'w') as h5f:
            for group in HT_format_dict:
                grp = h5f.create_group(group)
                for dataset in HT_format_dict[group]:
                    grp.create_dataset(dataset, data=HT_format_dict[group][dataset])



## Saving function

def monitor_to_h5(particles_to_save, filename, compression_opts=4,  readwrite_opts='a', verbose=False, collapse_into_avr=True):

        not_first_time = os.path.exists(filename)
        with h5py.File(filename, readwrite_opts) as hs:

            excluded_stuff=['XoStruct','extra_sources','add_particles']
            tracked_vars=['x','y','px','py']#,'delta','energy','kin_ps','kin_px','kin_py','kin_xprime','kin_yprime','ptau','pzeta','rpp','rvv','zeta'
            attrs=[a for a in dir(particles_to_save) if not (a.startswith('_')  or a in excluded_stuff )]

            for kk in attrs:

                if verbose: print('Writing '+kk)
                dat=getattr(particles_to_save , kk)
                
                if not_first_time:
                    if  kk in tracked_vars:
                        if collapse_into_avr:
                            dat=np.mean(dat,axis=0)
                        

                        hs[kk].resize( hs[kk].shape[0]+1, axis=0) 
                        hs[kk][-1,...] = dat

                else:
                    if isinstance(dat, np.ndarray ):

                        if collapse_into_avr:
                            dat=np.mean(dat,axis=0)        
                        maxshape=tuple([None]*len(np.array([dat]).shape))
                        
                        dset = hs.create_dataset(kk, shape=np.array([dat]).shape, dtype=dat.dtype, maxshape=maxshape , compression='gzip', compression_opts=compression_opts)
                        dset[...] = dat
                    else:
                        try:
                            hs[kk] = dat
                        except:
                            pass
        hs.close()




if __name__=="__main__":
    
    #study growth rate for negative chromaticity QPv 
    
    qpv_vect = np.linspace(0,-2,20)

    intensity_vect=[2.5e10] #np.linspace(3e10,8e10,6)
    #n_macroparticles / n_slices should be > 2000
    n_slices_0 = 300 #de 50 parriba (hasta 200)
    n_mp_0 = 1e6 #era 3e5

    storage_path='./forbench'


    for qpv in qpv_vect:
        for intensity in intensity_vect:
                print('--------------------------------------------------')
                print(f"Running chromaticity sweep, qpv : {qpv}...")
                run(qpv_val=qpv, n_slices=n_slices_0, n_mp=n_mp_0, intensity=intensity)
                
                if not os.path.exists(storage_path):
                    os.makedirs(storage_path)
                for f in glob.glob('./*.h5'):
                    shutil.move(f, storage_path)













import matplotlib.pylab as pl
import numpy as np
import os
import sys

#seed for random generation
pl.seed(0)

def push_simulation_to_folder(save_to_folder, data_from_folder):
    print "Copying all simulation results from %s " %data_from_folder\
          +"and simulation file to folder %s." % save_to_folder
    try:
        os.mkdir(save_to_folder)
    except(OSError):
        print "Result folder already exists. Overwriting..."
    os.system('cp %s/* %s' %(data_from_folder, save_to_folder))


def find_amplitude_at_comp(signal_with_spike, debug=False):
    n_compartments = len(signal_with_spike[:,0])
    amps = np.zeros(n_compartments)
    for comp in xrange(n_compartments):
        amps[comp] = np.max(signal_with_spike[comp,:])\
                     -np.min(signal_with_spike[comp,:])
        if debug:
            print amps[comp], ' [nA]'
    if debug:
        print "Number of compartments: %d" %n_compartments
        print "Max amplitude: %g" %np.max(amps)
    return amps

def find_comp_dist_from_elec(elec, cell, Mea):
    elec_x = Mea.elec_x
    elec_y = Mea.elec_y[elec]
    elec_z = Mea.elec_z[elec]
    n_compartments = len(cell.imem[:,0])
    dist = np.zeros(n_compartments)
    for comp in xrange(n_compartments):
        dist[comp] = np.sqrt((elec_x - cell.xmid[comp]/1000)**2 +\
                             (elec_y - cell.ymid[comp]/1000)**2 +\
                             (elec_z - cell.zmid[comp]/1000)**2)
    return dist

def find_comps_dist_from_soma(cell):
    n_compartments = len(cell.zmid)
    dist_from_soma = np.zeros(n_compartments)
    for comp in xrange(n_compartments):
        dist_from_soma[comp] = np.sqrt((cell.xmid[comp] - cell.xmid[0])**2\
                               + (cell.ymid[comp] - cell.ymid[0])**2\
                               + (cell.zmid[comp] - cell.zmid[0])**2)
    return dist_from_soma/1000

def single_comp_impact_on_elec(elec, mapping, cell, Mea):
    n_compartments = len(cell.xmid)
    for comp in xrange(n_compartments):
        impact_from_single_comp = mapping[elec, comp]
    return 0
   
def analyze_neuron(cell, mapping, Mea, signal_range):
    start_t, stop_t = signal_range
    start_t_ixd = np.argmin(np.abs(cell.tvec - start_t))
    stop_t_ixd  = np.argmin(np.abs(cell.tvec - stop_t))
    t_array = cell.tvec[start_t_ixd:stop_t_ixd]
    imem = cell.imem[:, start_t_ixd:stop_t_ixd]
    amps_at_comps = find_amplitude_at_comp(imem, debug=False)
    comp_dist_from_soma = find_comps_dist_from_soma(cell)
    #pl.figure()
    #pl.plot(comp_dist_from_soma, amps_at_comps, 'o')
    #pl.show()


    for elec in xrange(Mea.n_elecs):
        #if not elec == 110:
        #    continue
        print elec
        comp_dist_from_elec = find_comp_dist_from_elec(elec, cell, Mea)
        comp_impact_on_elec = mapping[elec,:]*amps_at_comps
        if 0:
            print np.argmax(comp_dist_from_elec)
            print np.argmin(comp_dist_from_elec)
            pl.figure()
            pl.axis('equal')
            pl.plot(cell.zmid/1000, cell.ymid/1000, 'o')
            pl.plot(Mea.elec_z, Mea.elec_y, 'o')
            pl.plot(cell.zmid[381]/1000, cell.ymid[381]/1000, 'D')
            pl.plot(cell.zmid[791]/1000, cell.ymid[791]/1000, 'D')
            pl.plot(Mea.elec_z[elec], Mea.elec_y[elec], 'x')
            pl.show() 
        pl.plot(comp_dist_from_elec, comp_impact_on_elec, 'o')
        
        pl.yscale('log')
        pl.xscale('log')
        pl.title('Electrode at distal tuft dendrite, when spike originates in tuft,\n'\
             +'and does not propagate to soma.')
        pl.xlabel('Compartment distance from electrode [mm]')
        pl.ylabel('Compartment impact on electrode [uV]')
    pl.show()
    
    
    #print amps_at_comps

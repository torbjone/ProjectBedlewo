import matplotlib.pylab as pl
import LFPy
import numpy as np
import os
import sys
import neuron
from plotting import plotstuff, simple_plot_2D,\
     plot_cell_compartments, plot_sec
from tools import push_simulation_to_folder, analyze_neuron, find_sec_from_comp
sim_folder = 'hay_model/'
LFPy.cell.neuron.load_mechanisms(sim_folder + '/mod')
pl.seed(999)
cellParameters = {
    'morphology' : sim_folder+'lfpy_version/morphologies/cell1.hoc',
    #'rm' : 30000,               # membrane resistance
    #'cm' : 1.0,                 # membrane capacitance
    'Ra' : 80,                 # axial resistance
    'v_init' : -80,             # initial crossmembrane potential
    #'e_pas' : -90,              # reversal potential passive mechs
    #'passive' : True,           # switch on passive mechs
    #'nsegs_method' : 'lambda_f',# method for setting number of segments,
    #'lambda_f' : 100,           # segments are isopotential at this frequency
    'timeres_NEURON' : 2**-4,   # dt of LFP and NEURON simulation.
    'timeres_python' : 2**-4,
    'tstartms' : -50,          #start time, recorders start at t=0
    'tstopms' : 120,           #stop time of simulation
    'custom_code'  : [sim_folder+'lfpy_version/custom_codes.hoc', \
                      sim_folder+'lfpy_version/biophys3.hoc'],
    # will if given list of files run this file
}

# Synaptic parameters taken from Hendrickson et al 2011
# Excitatory synapse parameters:
synapseParameters_AMPA = {
    'e' : 0,                    #reversal potential
    'syntype' : 'Exp2Syn',      #conductance based exponential synapse
    'tau1' : 1.,                #Time constant, rise
    'tau2' : 3.,                #Time constant, decay
    'weight' : 0.005,           #Synaptic weight
    'color' : 'r',              #for pl.plot
    'marker' : '.',             #for pl.plot
    'record_current' : True,    #record synaptic currents
}
# Excitatory synapse parameters
synapseParameters_NMDA = {         
    'e' : 0,
    'syntype' : 'Exp2Syn',
    'tau1' : 10.,
    'tau2' : 30.,
    'weight' : 0.005,
    'color' : 'm',
    'marker' : '.',
    'record_current' : True,
}
# Inhibitory synapse parameters
synapseParameters_GABA_A = {         
    'e' : -80,
    'syntype' : 'Exp2Syn',
    'tau1' : 1.,
    'tau2' : 12.,
    'weight' : 0.005,
    'color' : 'b',
    'marker' : '.',
    'record_current' : True
}
# where to insert, how many, and which input statistics
insert_synapses_AMPA_args = {
    'section' : 'apic',
    'n' : 100,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellParameters['tstartms'], cellParameters['tstopms'], 2, 10]
}
insert_synapses_NMDA_args = {
    'section' : 'alldend',
    'n' : 40,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellParameters['tstartms'], cellParameters['tstopms'], 5, 20]
}
insert_synapses_GABA_A_args = {
    'section' : 'dend',
    'n' : 100,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellParameters['tstartms'], cellParameters['tstopms'], 2, 10]
}

clamp_1 = {
    'idx' : 0,
    'record_current' : True,
    'amp' : 1.9, #[nA]
    'dur' : 5.,
    'delay' :20,
    #'freq' : 10,
    #'phase' : 0,
    #'pkamp' : 300e-3,
    'pptype' : 'IClamp',
}
clamp_2 = {
    'idx' : 611,
    'record_current' : True,
    'amp' : 1.9, #[nA]
    'dur' : 1,
    'delay' :23,
    #'freq' : 10,
    #'phase' : 0,
    #'pkamp' : 300e-3,
    'pptype' : 'IClamp',
}
# Parameters for the cell.simulate() call, recording membrane- and syn.-currents
simulationParameters = {
    'rec_imem' : True,  # Record Membrane currents during simulation
    'rec_isyn' : True,  # Record synaptic currents
    'rec_vmem' : True,    #record membrane potential for all compartments
    #'rec_variables': ['ina_NaTa_t'],
}

def get_cell(output_folder, do_simulation = True):
    def insert_synapses(synparams, section, n, spTimesFun, args):
        #find n compartments to insert synapses onto
        #idx = [clamp_1['idx'], clamp_2['idx']]#
        idx = cell.get_rand_idx_area_norm(section=section, nidx=n)
        #spiketimes = [np.array([15, 25]), np.array([16,26])]
        #Insert synapses in an iterative fashion
        for count, index in enumerate(idx):
            synparams.update({'idx' : int(index)})
            s = LFPy.Synapse(cell,**synparams)
            spiketimes = spTimesFun(args[0], args[1], args[2], args[3])
            s.set_spike_times(spiketimes)
            #s.set_spike_times(spiketimes[count])
    cell = LFPy.Cell(**cellParameters)
    cell.set_rotation(x = -pl.pi/2)
    cell.set_rotation(y = pl.pi/2)
    #find_sec_from_comp(cell, 756)
    cell.set_pos(xpos = -25)
    #cell.set_rotation(x = pl.pi/2)
    #cell.set_rotation(y = -pl.pi/2)
    #cell.set_rotation(y = pl.pi/10)
    cell.set_rotation(z = -pl.pi/165)


    
    #find_sec_from_comp(cell, 644)
    #find_sec_from_comp(cell, 801)
    #plot_sec(cell, 'apic[76]')
    #sys.exit()
    def insert_glutamate_stim(cell, section = 'apic[63]', site = 0.5):
        gmaxS=10
        neuron.h('access %s' %section)
        glut_syn = neuron.h.glutamate(site)
        glut_syn.delay = 50
        glut_syn.ntar = 1.3
        glut_syn.gmax = gmaxS
        glut_syn.Nspike=3
        glut_syn.Tspike=20
        return glut_syn

    def insert_many_glutamate_stim(cell, nSyn = 500):
	import random
        gmaxS=1
	syn_array = []
	i = 0
	while i < nSyn:
	    apic_sec = random.randint(0,108)
	    section = 'apic[%d]' % apic_sec
	    comps = cell.get_idx_section(section)
	    idx = random.choice(comps)
	    if cell.ymid[idx] > 600:
	        #syns_positions.append(
		neuron.h('access %s' % section)
		syn_array.append(neuron.h.glutamate(random.uniform(0.1,0.95)))
		syn_array[-1].delay = 50
		syn_array[-1].ntar = 1
		syn_array[-1].gmax = gmaxS
		syn_array[-1].Nspike=3
		syn_array[-1].Tspike=20
		i += 1
        return syn_array#, syns_positions
    
    if do_simulation:
        os.system('cp %s %s' %(sys.argv[0], output_folder))
        np.save(output_folder + 'x_start.npy', cell.xstart)
        np.save(output_folder + 'y_start.npy', cell.ystart)
        np.save(output_folder + 'z_start.npy', cell.zstart)
        np.save(output_folder + 'x_end.npy', cell.xend)
        np.save(output_folder + 'y_end.npy', cell.yend)
        np.save(output_folder + 'z_end.npy', cell.zend)
        np.save(output_folder + 'diam.npy', cell.diam)

        syn1 = insert_glutamate_stim(cell, \
                section = find_sec_from_comp(cell, 615), site = 0.5)
        syn2 = insert_glutamate_stim(cell, \
                section = find_sec_from_comp(cell, 646), site = 0.59)
        syn3 = insert_glutamate_stim(cell, \
                section = find_sec_from_comp(cell, 671) , site = 0.5)
        syn4 = insert_glutamate_stim(cell, \
                 section = find_sec_from_comp(cell, 623), site = 0.2)
        syn5 = insert_glutamate_stim(cell, \
                 section = find_sec_from_comp(cell, 628), site = 0.9)
        syn6 = insert_glutamate_stim(cell, \
                section = find_sec_from_comp(cell, 657), site = 0.8) 
        #currentClamp_1 = LFPy.StimIntElectrode(cell, **clamp_1)
        #currentClamp_2 = LFPy.StimIntElectrode(cell, **clamp_2)
        #insert_synapses(synapseParameters_AMPA, **insert_synapses_AMPA_args)
        #insert_synapses(synapseParameters_NMDA, **insert_synapses_NMDA_args)
        #insert_synapses(synapseParameters_GABA_A, **insert_synapses_GABA_A_args)
        #syn_array = insert_many_glutamate_stim(cell, nSyn = 100)
        #neuron.h('objref n_vec')
        #neuron.h('n_vec = new Vector()')
        #neuron.h('n_vec.record(&soma.m_hh(0.5))')
        
        cell.simulate(**simulationParameters)
        np.save(output_folder + 'imem.npy', cell.imem)
        np.save(output_folder + 'vmem.npy', cell.vmem)
        np.save(output_folder + 'tvec.npy', cell.tvec)
        plotstuff(cell, clamp_1, clamp_2)
        os.system('cp %s %s' %('example_fig.png', output_folder))
    else:
        cell.vmem = np.load(output_folder + 'vmem.npy')
        cell.imem = np.load(output_folder + 'imem.npy')
        cell.tvec = np.load(output_folder + 'tvec.npy')
    return cell

if __name__ == '__main__':
    output_folder = 'hay_results/initial_test/'
    do_simulation = True
    plot_range = [49,140]
    try:
        os.mkdir(output_folder)
    except(OSError):
        if do_simulation:
            print "Result folder already exists. Overwriting..."
        else:
            print "Loading simulation files..."
    cell = get_cell(output_folder, do_simulation)
    
    #plot_cell_compartments(cell)
    #simple_plot_2D(cell, plot_range, clamp_1, clamp_2)
    #push_simulation_to_folder('extracellular_test/', 'hay_results/initial_test/')
   

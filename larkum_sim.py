import matplotlib.pylab as pl
import LFPy
import numpy as np
import os
import sys
import neuron
from plotting import plotstuff, simple_plot_2D,\
     plot_cell_compartments
sim_folder = 'larkum_model/'
#load compiled mechs from the mod-folder
LFPy.cell.neuron.load_mechanisms(sim_folder + '/mod')


cellParameters = {
    'morphology' : '070603c2_copy.hoc',
    'rm' : 30000,               # membrane resistance
    'cm' : 1.0,                 # membrane capacitance
    'Ra' : 150,                 # axial resistance
    'v_init' : -65,             # initial crossmembrane potential
    'e_pas' : -65,              # reversal potential passive mechs
    'passive' : True,           # switch on passive mechs
    'nsegs_method' : 'lambda_f',# method for setting number of segments,
    'lambda_f' : 100,           # segments are isopotential at this frequency
    'timeres_NEURON' : 2**-5,   # dt of LFP and NEURON simulation.
    'timeres_python' : 2**-5,
    'tstartms' : -100,          #start time, recorders start at t=0
    'tstopms' : 200,           #stop time of simulation
    'custom_code'  : ['apical_simulation.hoc'],        # will if given list of files run this file
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
    'n' : 20,
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
    'idx' : 202,
    'record_current' : True,
    'amp' : 10.0, #[nA]
    'dur' : 0.5,
    'delay' :150,
    #'freq' : 10,
    #'phase' : 0,
    #'pkamp' : 300e-3,
    'pptype' : 'IClamp',
}
clamp_2 = {
    'idx' : 420,
    'record_current' : True,
    'amp' : 1.3, #[nA]
    'dur' : 1,
    'delay' :27,
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
}

def get_cell(output_folder, do_simulation = True):
    def insert_synapses(synparams, section, n, spTimesFun, args):
        #find n compartments to insert synapses onto
        #idx = [clamp_1['idx'], clamp_2['idx']]#
        idx = cell.get_rand_idx_area_norm(section=section, nidx=n)

        #spiketimes = [np.array([15, 25]), np.array([16,26])]
        #Insert synapses in an iterative fashion
        for count, i in enumerate(idx):
            synparams.update({'idx' : int(i)})
            s = LFPy.Synapse(cell,**synparams)
            spiketimes = spTimesFun(args[0], args[1], args[2], args[3])
            s.set_spike_times(spiketimes)
            #s.set_spike_times(spiketimes[count])
    cell = LFPy.Cell(**cellParameters)
    cell.set_rotation(x = pl.pi/2)
    #cell.set_rotation(z = pl.pi/2)
    if do_simulation:
        os.system('cp %s %s' %(sys.argv[0], output_folder))
        np.save(output_folder + 'x_start.npy', cell.xstart)
        np.save(output_folder + 'y_start.npy', cell.ystart)
        np.save(output_folder + 'z_start.npy', cell.zstart)
        np.save(output_folder + 'x_end.npy', cell.xend)
        np.save(output_folder + 'y_end.npy', cell.yend)
        np.save(output_folder + 'z_end.npy', cell.zend)
        np.save(output_folder + 'diam.npy', cell.diam)
        currentClamp_1 = LFPy.StimIntElectrode(cell, **clamp_1)
        #currentClamp_2 = LFPy.StimIntElectrode(cell, **clamp_2)
        insert_synapses(synapseParameters_AMPA, **insert_synapses_AMPA_args)
        insert_synapses(synapseParameters_NMDA, **insert_synapses_NMDA_args)
        insert_synapses(synapseParameters_GABA_A, **insert_synapses_GABA_A_args)
        cell.simulate(**simulationParameters)
        np.save(output_folder + 'imem.npy', cell.imem)
        np.save(output_folder + 'vmem.npy', cell.vmem)
        np.save(output_folder + 'tvec.npy', cell.tvec)
        plotstuff(cell)
        os.system('cp %s %s' %('example_fig.png', output_folder))
    else:
        cell.vmem = np.load(output_folder + 'vmem.npy')
        cell.imem = np.load(output_folder + 'imem.npy')
        cell.tvec = np.load(output_folder + 'tvec.npy')
    return cell


if __name__ == '__main__':
    output_folder = 'larkum_results/initial_test/'
    do_simulation = True
    plot_range = [149,155]
    try:
        os.mkdir(output_folder)
    except(OSError):
        if do_simulation:
            print "Result folder already exists. Overwriting..."
        else:
            print "Loading simulation files..."
    cell = get_cell(output_folder, do_simulation)
    
    plot_cell_compartments(cell)
    simple_plot_2D(cell, plot_range, clamp_1, clamp_2)
    #push_simulation_to_folder('two_syn_AP_in_soma_long_delay', output_folder)

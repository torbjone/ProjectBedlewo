import matplotlib.pylab as pl
import LFPy
import numpy as np
import os
import sys
import neuron
from plotting import plotstuff, simple_plot_2D,\
     plot_cell_compartments
from tools import push_simulation_to_folder, analyze_neuron

class HOCProxy(object):
	def __init__(self, interp):
		self.interp = interp

	def __getattr__(self, name,):
		fnc = getattr(self.interp, name)
		def my_call(*args,**kwargs):
			print name, args, kwargs
			return fnc(*args,**kwargs)
		return my_call
			

	def __call__(self, *args, **kwargs):
		print args, kwargs
		return self.interp(*args, **kwargs)

old_interp = neuron.h
#neuron.h = HOCProxy(old_interp)


sim_folder = 'larkum_model/'
LFPy.cell.neuron.load_mechanisms(sim_folder)

cellParameters = {
    'morphology' : sim_folder + '070603c2_copy.hoc',
    #'rm' : 20000,               # membrane resistance
    #'cm' : 1.,                 # membrane capacitance
    'Ra' : 80,                 # axial resistance
    #'v_init' : -65,             # initial crossmembrane potential
    #'e_pas' : -70,              # reversal potential passive mechs
    'passive' : False,           # switch on passive mechs
    #'nsegs_method' : 'lambda_f',# method for setting number of segments,
    #'lambda_f' : 500,           # segments are isopotential at this frequency
    'timeres_NEURON' : 0.125/2,   # dt of LFP and NEURON simulation.
    'timeres_python' : 0.125/2,
    'tstartms' : 0,          #start time, recorders start at t=0
    'tstopms' : 120,           #stop time of simulation
    'custom_code'  : [ sim_folder +'apical_simulation_changed.hoc'],     # will if given list of files run this file
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
    'name' : 'AMPA',
    'section' : 'allsec',
    'n' : 100,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellParameters['tstartms'], cellParameters['tstopms'], 2, 10]
}
insert_synapses_NMDA_args = {
    'name' : 'NMDA',
    'section' : 'alldend',
    'n' : 20,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellParameters['tstartms'], cellParameters['tstopms'], 5, 20]
}
insert_synapses_GABA_A_args = {
    'name' : 'GABA_A',
    'section' : 'dend',
    'n' : 100,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellParameters['tstartms'], cellParameters['tstopms'], 2, 10]
}

clamp_1 = {
    'idx' : 142,
    'record_current' : True,
    'amp' : 0.6, #[nA]
    'dur' : 400,
    'delay' :50,
    'pptype' : 'IClamp',
}

clamp_2 = {
    'idx' : 349,
    'record_current' : True,
    'amp' : 0.8, #[nA]
    'dur' : 2,
    'delay' :21,
    'pptype' : 'IClamp',
}

clamp_3 = {
    'idx' : 0,
    'record_current' : True,
    'amp' : .3, #[nA]
    'dur' : 30,
    'delay' :10,
    'pptype' : 'IClamp',
}


# Parameters for the cell.simulate() call, recording membrane- and syn.-currents
simulationParameters = {
    'rec_imem' : True,  # Record Membrane currents during simulation
    'rec_isyn' : True,  # Record synaptic currents
    'rec_vmem' : True,    #record membrane potential for all compartments
    'rec_variables' : ['ica', 'cai']
}

def get_cell(output_folder, do_simulation = True):
    def insert_synapses(synparams, name, section, n, spTimesFun, args):
        #find n compartments to insert synapses onto
        #idx = [clamp_1['idx'], clamp_2['idx']]#
        idx = cell.get_rand_idx_area_norm(section=section, nidx=n)
        #if name == 'NMDA':
        #    idx_tuft = []
       #     for i in idx:
       #         if cell.zmid[i] > 500:
       #             idx_tuft.append(i)
       #     idx = idx_tuft
        #spiketimes = [np.array([15, 25]), np.array([16,26])]
        #Insert synapses in an iterative fashion
        for count, i in enumerate(idx):
            synparams.update({'idx' : int(i)})
            s = LFPy.Synapse(cell,**synparams)
            spiketimes = spTimesFun(args[0], args[1], args[2], args[3])
            s.set_spike_times(spiketimes)
            #s.set_spike_times(spiketimes[count])

    def insert_glutamate_stim(cell, section = 'apic[63]', site = 0.5):
        gmaxS=20
        neuron.h('access %s' %section)
        glut_syn = neuron.h.glutamate(site)
        glut_syn.delay = 50
        glut_syn.ntar = 1
        glut_syn.gmax = gmaxS
        glut_syn.Nspike=2
        glut_syn.Tspike=20
        return glut_syn

    def insert_many_glutamate_stim(cell, nSyn = 500):
	import random
        gmaxS=1
	syn_array = []
	i = 0
	while i < nSyn:
	    apic_sec = random.randint(0,115)
	    section = 'apic[%d]' % apic_sec
	    comps = cell.get_idx_section(section)
	    idx = random.choice(comps)
	    if cell.ymid[idx] > 800:
	        #syns_positions.append(
		neuron.h('access %s' % section)
		syn_array.append(neuron.h.glutamate(random.uniform(0.5,0.95)))
		syn_array[-1].delay = 50
		syn_array[-1].ntar = 1
		syn_array[-1].gmax = gmaxS
		syn_array[-1].Nspike=3
		syn_array[-1].Tspike=20
		i += 1
        return syn_array#, syns_positions

    cell = LFPy.Cell(**cellParameters)
    cell.apic_11 = cell.get_closest_idx(x=52.52, y=441.92, z=0.75)
    cell.apic_58 = cell.get_closest_idx(x=29.05, y=742.84, z=-34.05)
    cell.apic_59 = cell.get_closest_idx(x = 20.9, y = 825.51, z = -40.25)
  
    cell.set_pos(xpos = -45)
    cell.set_rotation(y = -pl.pi/2)
    cell.set_rotation(y = pl.pi/10)
    cell.set_rotation(z = +pl.pi/65)
    
    if do_simulation:
        os.system('cp %s %s' %(sys.argv[0], output_folder))
        np.save(output_folder + 'x_start.npy', cell.xstart)
        np.save(output_folder + 'y_start.npy', cell.ystart)
        np.save(output_folder + 'z_start.npy', cell.zstart)
        np.save(output_folder + 'x_end.npy', cell.xend)
        np.save(output_folder + 'y_end.npy', cell.yend)
        np.save(output_folder + 'z_end.npy', cell.zend)
        np.save(output_folder + 'diam.npy', cell.diam)
        try:
             np.save(output_folder + 'ica.npy',cell.rec_variables['ica'])
             np.save(output_folder + 'cai.npy',cell.rec_variables['cai'])
        except:
            pass
        #currentClamp_1 = LFPy.StimIntElectrode(cell, **clamp_1)
        #currentClamp_2 = LFPy.StimIntElectrode(cell, **clamp_2)
        #currentClamp_3 = LFPy.StimIntElectrode(cell, **clamp_3)
        #insert_synapses(synapseParameters_AMPA, **insert_synapses_AMPA_args)
        #insert_synapses(synapseParameters_NMDA, **insert_synapses_NMDA_args)
	#insert_synapses(synapseParameters_GABA_A, **insert_synapses_GABA_A_args)
        #glut_syn = insert_glutamate_stim(cell, section = 'apic[63]')
	#glut_syn2 = insert_glutamate_stim(cell, section = 'apic[45]')
	#glut_syn = insert_many_glutamate_stim(cell, 100)
        #glut_syn1 = insert_glutamate_stim(cell, section = 'apic[63]', site = 0.59)
        #glut_syn2 = insert_glutamate_stim(cell, section = 'apic[71]', site = 0.59)
	#glut_syn3 = insert_glutamate_stim(cell, section = 'apic[71]', site = 0.31)
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
        try:
             cell.rec_variables['ica'] = np.load(output_folder + 'ica.npy')
             cell.rec_variables['cai'] = np.load(output_folder + 'cai.npy')
        except:
            pass
    return cell

if __name__ == '__main__':
    #output_folder = 'larkum_results/initial_test/'
    output_folder = 'larkum_no_TTX_superthreshold/'
    #output_folder = 'larkum_Spikes_traveling_both_ways_initiated_by_inj/'
    do_simulation = True
    plot_range = [49,120]
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
    #push_simulation_to_folder('larkum_TTX_many_synapses/', output_folder)
    #analyze_neuron(cell, signal_range = [22,26])



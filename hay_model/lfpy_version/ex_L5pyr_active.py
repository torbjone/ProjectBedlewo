#!/usr/bin/env python

# importing some modules, setting some matplotlib values for pl.plot.
import pylab as pl
import LFPy

#load compiled mechs from the mod-folder
LFPy.cell.neuron.load_mechanisms("../mod")

pl.rcParams.update({'font.size' : 10, 'figure.figsize' : [16,9],'wspace' : 0.5 ,'hspace' : 0.5})

#seed for random generation
pl.seed(9876543210)

#plot pops up by itself
pl.interactive(1)

################################################################################
# A couple of function declarations
################################################################################

def plotstuff():
    fig = pl.figure(figsize=[12, 8])
    
    ax = fig.add_axes([0.1, 0.7, 0.5, 0.2])
    ax.plot(cell.tvec,cell.somav)
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('Soma pot. [mV]')
    
    ax = fig.add_axes([0.1, 0.4, 0.5, 0.2])
    for i in xrange(len(cell.synapses)):
        ax.plot(cell.tvec,cell.synapses[i].i,color=cell.synapses[i].color)
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('Syn. i [nA]')
    
    ax = fig.add_axes([0.1, 0.1, 0.5, 0.2])
    absmaxLFP = abs(pl.array([electrode.LFP.max(),electrode.LFP.min()])).max()
    pl.imshow(electrode.LFP,vmax=absmaxLFP,vmin=-absmaxLFP,origin='lower',
           extent=(cell.tvec[0],cell.tvec[-1],electrode.z[0],electrode.z[-1]),cmap='jet_r',
           interpolation='nearest')
    cbar = pl.colorbar(ax=ax)
    cbar.set_label('LFP (mV)')
    pl.axis('tight')
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('z [$\mu$m]')
    
    ax = fig.add_axes([0.65, 0.1, 0.25, 0.8], frameon=False)
    for i in xrange(cell.xend.size):
        ax.plot([cell.xstart[i],cell.xend[i]],[cell.zstart[i],cell.zend[i]],color='k')
    for i in xrange(len(cell.synapses)):
        ax.plot([cell.synapses[i].x],[cell.synapses[i].z],\
            color=cell.synapses[i].color,marker=cell.synapses[i].marker)
    for i in xrange(electrode.x.size):
        ax.plot(electrode.x[i],electrode.z[i],color='g',marker='o')
    pl.axis('equal')
    pl.axis(pl.array(pl.axis())*0.8)
    ax.set_xticks([])
    ax.set_yticks([])

def insert_synapses(cell, synparams, section, n, spTimesFun, args):
    #find n compartments to insert synapses onto
    idx = cell.get_rand_idx_area_norm(section=section, nidx=n)

    #Insert synapses in an iterative fashion
    for i in idx:
        synparams.update({'idx' : int(i)})

        # Some input spike train using the function call
        spiketimes = spTimesFun(args[0], args[1], args[2], args[3])

        # Create synapse(s) and setting times using the Synapse class in LFPy
        synapse = LFPy.Synapse(cell,**synparams)
        synapse.set_spike_times(spiketimes)

################################################################################
# Define parameters, using dictionaries
# It is possible to set a few more parameters for each class or functions, but
# we chose to show only the most important ones here.
################################################################################

cellparams = {          #various cell parameters,
    'morphology' : 'morphologies/cell1.asc', #Mainen&Sejnowski, Nature, 1996
    'v_init' : -80,     #initial crossmembrane potential
    'passive' : False,   #switch on passive mechs
    'nsegs_method' : 'lambda_f',
    'lambda_f' : 100,
    'timeres_NEURON' : 2**-4,   #[ms] dt's should be in powers of 2 for both,
    'timeres_python' : 2**-4,   #need binary representation
    'tstartms' : -200,  #start time of simulation, recorders start at t=0
    'tstopms' : 1000,   #stop simulation at 1000 ms. these can be overridden
                        #by setting these arguments in cell.simulation()
    'custom_code'  : ['custom_codes.hoc', 'biophys1.hoc'],    #Custom .hoc/.py-scripts
}
#Synaptic parameters taken from Hendrickson et al 2011
synparams_AMPA = {         #Excitatory synapse parameters
    'e' : 0,           #reversal potential
    'syntype' : 'Exp2Syn',   #conductance based exponential synapse
    'tau1' : 1.,         #Time constant, rise
    'tau2' : 3.,         #Time constant, decay
    'weight' : 0.005,   #Synaptic weight
    'color' : 'r',      #for pl.plot
    'marker' : '.',     #for pl.plot
    'record_current' : True,    #record synaptic currents
}
synparams_NMDA = {         #Excitatory synapse parameters
    'e' : 0,           #reversal potential
    'syntype' : 'Exp2Syn',   #conductance based exponential synapse
    'tau1' : 10.,         #Time constant, rise
    'tau2' : 30.,         #Time constant, decay
    'weight' : 0.005,   #Synaptic weight
    'color' : 'm',      #for pl.plot
    'marker' : '.',     #for pl.plot
    'record_current' : True,    #record synaptic currents
}
synparams_GABA_A = {         #Inhibitory synapse parameters
    'e' : -80,
    'syntype' : 'Exp2Syn',
    'tau1' : 1.,
    'tau2' : 12.,
    'weight' : 0.005,
    'color' : 'b',
    'marker' : '.',
    'record_current' : True
}
#where to insert, how many, and which input statistics
insert_synapses_AMPA_args = {
    'section' : 'apic',
    'n' : 200,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellparams['tstartms'], cellparams['tstopms'], 2, 10]
}
insert_synapses_NMDA_args = {
    'section' : 'alldend',
    'n' : 10,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellparams['tstartms'], cellparams['tstopms'], 5, 20]
}
insert_synapses_GABA_A_args = {
    'section' : 'dend',
    'n' : 100,
    'spTimesFun' : LFPy.inputgenerators.stationary_gamma,
    'args' : [cellparams['tstartms'], cellparams['tstopms'], 2, 10]
}

N = pl.empty((16, 3))
for i in xrange(N.shape[0]): N[i,] = [1, 0, 0] #normal unit vec. to contacts
electrodeparams = {             #parameters for electrode class
    'sigma' : 0.3,              #Extracellular potential
    'x' : pl.zeros(16)+10,      #Coordinates of electrode contacts
    'y' : pl.zeros(16),
    'z' : pl.linspace(-500,1000,16),
    'n' : 20,
    'r' : 10,
    'N' : N,
}


################################################################################
# Main simulation procedure
################################################################################

#pl.close('all')   #close open figures

#Initialize cell instance, using the LFPy.Cell class
cell = LFPy.Cell(**cellparams)

#Insert synapses
insert_synapses(cell, synparams_AMPA, **insert_synapses_AMPA_args)
insert_synapses(cell, synparams_NMDA, **insert_synapses_NMDA_args)
insert_synapses(cell, synparams_GABA_A, **insert_synapses_GABA_A_args)

electrode = LFPy.RecExtElectrode(**electrodeparams)
#perform NEURON simulation, results saved as attributes in the cell instance
cell.simulate(electrode, rec_isyn = True)

#plotting some variables and geometry, saving output to .pdf.
plotstuff()
################################################################################

#!/usr/bin/env python
################################################################################
#
# This is an example scripts using LFPy with an active cell model adapted from
# Mainen and Sejnowski, Nature 1996, for the original files, see
# http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=2488
#
# This scripts is set up to use the model, where the active conductances are set
# in the file "active_declarations_example2.hoc", and uses the mechanisms from
# the .mod-files provided here. For this example to work, run "nrnivmodl" in
# this folder to compile these mechanisms
# (i.e. /$PATHTONEURON/nrn/x86_64/bin/nrnivmodl).
#
# A single excitatory synapse drive the neuron into producing a single action-
# potential, and the local field potential are calculated on a dense 2D-grid
# on around the soma.
#
################################################################################
#import some plotting stuff and the LFPy-module
import matplotlib.pylab as pl
from mayavi import mlab
import LFPy
import neuron
from time import time
import ipdb

#load compiled mechs from the mod-folder
LFPy.cell.neuron.load_mechanisms("../mod")

#set some plotting parameters
pl.rcParams.update({'font.size' : 15,
    'figure.facecolor' : '1',
    'left': 0.1, 'wspace' : 0.5 ,'hspace' : 0.5})


################################################################################
# Function declarations
################################################################################

def plotstuff(cell, electrode):

    figure = mlab.figure(size=(800,800))
    
    l_list = []
    for sec in neuron.h.allsec():
        idx = cell.get_idx_section(sec.name())
        j = 0
        for seg in sec:
            i = idx[j]
            x = pl.array([cell.xstart[i],cell.xend[i]])
            y = pl.array([cell.ystart[i],cell.yend[i]])
            z = pl.array([cell.zstart[i],cell.zend[i]])
            s = pl.array([seg.v, seg.v])
            
            l = mlab.plot3d(x, y, z, s, colormap = 'Spectral',
                            tube_radius=cell.diam[i],
                            representation='surface', vmin=-70, vmax=10)
            l_list.append(l)
            print j
            j += 1
    
    t0 = time()
    ipdb.set_trace()    
    #ms = l_list[0].mlab_source
    while time()-t0 < 10:
        for l in l_list:
            ms = l.mlab_source
            s = pl.rand()*80 -70
            scalars = pl.array([s, s])
            ms.set(scalars = scalars)
        
    
    
    
    #
    #for i in xrange(cell.totnsegs):
    #    x = pl.array([cell.xstart[i],cell.xend[i]])
    #    y = pl.array([cell.ystart[i],cell.yend[i]])
    #    z = pl.array([cell.zstart[i],cell.zend[i]])
    #    
    #    l = mlab.plot3d(x, y, z, tube_radius=cell.diam[i],
    #                    color=(1,1,1), representation='surface')
    #    #figure=gcf())
    #p = mlab.points3d(cell.xstart, cell.ystart, cell.zstart, cell.diam, color=(1,1,1), scale_factor=2.0)
    #p = mlab.points3d(cell.xend, cell.yend, cell.zend, cell.diam, color=(1,1,1), scale_factor=2.0)
    #
    
    
    
    
    #
    ##creating array of points and corresponding diameters along structure
    #xcoords = pl.array([])
    #ycoords = pl.array([])
    #zcoords = pl.array([])
    #diams = pl.array([])
    #    
    #for i in xrange(cell.xend.size):
    #    xcoords = pl.r_[xcoords, pl.linspace(cell.xstart[i],
    #                                cell.xend[i], int(cell.length[i]/cell.diam[i]))]   
    #    ycoords = pl.r_[ycoords, pl.linspace(cell.ystart[i],
    #                                cell.yend[i], int(cell.length[i]/cell.diam[i]))]   
    #    zcoords = pl.r_[zcoords, pl.linspace(cell.zstart[i],
    #                                cell.zend[i], int(cell.length[i]/cell.diam[i]))]   
    #    diams = pl.r_[diams, pl.linspace(cell.diam[i], cell.diam[i],
    #                                int(cell.length[i]/cell.diam[i]))]
    #
    ##sort along depth-axis
    #argsort = pl.argsort(ycoords)
    #
    #
    #
    #
    #mlab.points3d(xcoords[argsort], ycoords[argsort], zcoords[argsort], diams[argsort],
    #              colormap="copper", scale_factor=2.0)
    #
    
    ##plotting
    #fig = pl.figure(figsize=[15, 10])
    #ax = fig.add_axes([0.1, 0.1, 0.533334, 0.8], frameon=False)
    #ax.scatter(xcoords[argsort], zcoords[argsort], s=diams[argsort]**2*20,
    #           c=ycoords[argsort], edgecolors='none', cmap='gray')
    #ax.plot(electrode.x, electrode.z, '.', marker='o', markersize=5, color='k')
    #
    #i = 0
    #limLFP = abs(electrode.LFP).max()
    #for LFP in electrode.LFP:
    #    tvec = cell.tvec*0.6 + electrode.x[i] + 2
    #    if abs(LFP).max() >= 1:
    #        factor = 2
    #        color='r'
    #    elif abs(LFP).max() < 0.25:
    #        factor = 50
    #        color='b'
    #    else:
    #        factor = 10
    #        color='g'
    #    trace = LFP*factor + electrode.z[i]
    #    ax.plot(tvec, trace, color=color, lw = 2)
    #    i += 1
    #
    #ax.plot([22, 28], [-60, -60], color='k', lw = 3)
    #ax.text(22, -65, '10 ms')
    #
    #ax.plot([40, 50], [-60, -60], color='k', lw = 3)
    #ax.text(42, -65, '10 $\mu$m')
    #
    #ax.plot([60, 60], [20, 30], color='r', lw=2)
    #ax.text(62, 20, '5 mV')
    #
    #ax.plot([60, 60], [0, 10], color='g', lw=2)
    #ax.text(62, 0, '1 mV')
    #
    #ax.plot([60, 60], [-20, -10], color='b', lw=2)
    #ax.text(62, -20, '0.1 mV')
    #
    #
    #
    #ax.set_xticks([])
    #ax.set_yticks([])
    #
    #ax.axis([-61, 61, -61, 61])
    #
    #ax.set_title('Location-dependent extracellular spike shapes')
    #
    ##plotting the soma trace    
    #ax = fig.add_axes([0.75, 0.55, 0.2, 0.35])
    #ax.plot(cell.tvec, cell.somav)
    #ax.set_title('Somatic action-potential')
    #ax.set_ylabel(r'$V_\mathrm{membrane}$ (mV)')
    #
    ##plotting the synaptic current
    #ax = fig.add_axes([0.75, 0.1, 0.2, 0.35])
    #ax.plot(cell.tvec, cell.synapses[0].i)
    #ax.set_title('Synaptic current')
    #ax.set_ylabel(r'$i_\mathrm{synapse}$ (nA)')
    #ax.set_xlabel(r'time (ms)')

################################################################################
# Define parameters, using dictionaries
# It is possible to set a few more parameters for each class or functions, but
# we chose to show only the most important ones here.
################################################################################

#define cell parameters used as input to cell-class
cellParameters = {
    'morphology' : 'morphologies/cell2.hoc',
    #'rm' : 30000,               # membrane resistance
    #'cm' : 1.0,                 # membrane capacitance
    #'Ra' : 150,                 # axial resistance
    #'v_init' : -65,             # initial crossmembrane potential
    #'e_pas' : -65,              # reversal potential passive mechs
    'passive' : False,           # switch on passive mechs
    'nsegs_method' : 'lambda_f',# method for setting number of segments,
    'lambda_f' : 1,           # segments are isopotential at this frequency
    'timeres_NEURON' : 2**-5,   # dt of LFP and NEURON simulation.
    'timeres_python' : 2**-5,
    'tstartms' : -100,           #start time, recorders start at t=0
    'tstopms' : 10,             #stop time of simulation
    'custom_code' : ['custom_codes.hoc', 'biophys1.hoc'],
}

#Synaptic parameters, corresponding to a NetCon synapse built into NEURON
synapseParameters = {
    'idx' : 0,               # insert synapse on index "0", the soma
    'e' : 0.,                # reversal potential of synapse
    'syntype' : 'Exp2Syn',   # conductance based double-exponential synapse
    'tau1' : 1.0,            # Time constant, rise
    'tau2' : 1.0,            # Time constant, decay
    'weight' : 0.05,         # Synaptic weight
    'record_current' : True, # Will enable synapse current recording
}

#Generate the grid in xz-plane over which we calculate local field potentials
x = pl.linspace(-50, 50, 11)
z = pl.linspace(-50, 50, 11)
X, Z = pl.meshgrid(x, z)
y = pl.zeros(X.size)

#define parameters for extracellular recording electrode, using optional method
electrodeParameters = {
    'sigma' : 0.3,              # extracellular conductivity
    'x' : X.reshape(-1),        # x,y,z-coordinates of contact points
    'y' : y,
    'z' : Z.reshape(-1),
    'method' : 'som_as_point',  #treat soma segment as sphere source
}

################################################################################
# Main simulation procedure, setting up extracellular electrode, cell, synapse
################################################################################

#close open figures
pl.close('all')

#create extracellular electrode object
electrode = LFPy.RecExtElectrode(**electrodeParameters)

#Initialize cell instance, using the LFPy.Cell class
cell = LFPy.Cell(**cellParameters)
#set the position of midpoint in soma to Origo (not needed, this is the default)
cell.set_pos(xpos = 0, ypos = 0, zpos = 0)
#rotate the morphology 90 degrees around z-axis
cell.set_rotation(z = pl.pi/2)

#attach synapse with parameters and set spike time
synapse = LFPy.Synapse(cell, **synapseParameters)
synapse.set_spike_times(pl.array([1]))

#perform NEURON simulation, results saved as attributes in the cell instance
#cell.simulate(electrode = electrode, rec_isyn=True)

# Plotting of simulation results:
plotstuff(cell, electrode)
#pl.savefig('example2.pdf')

#pl.show()



import pylab as pl
import numpy as np

def plot_interval_on_all_elecs(cell, signal_at_elec, t_interval, Mea):
    pl.close('all')
    n_elecs = Mea.n_elecs
    n_compartments = len(cell.imem[:,0])
    electrode_separation = Mea.electrode_separation
    t_array = cell.tvec
    t_start, t_stop = t_interval[:]
    t_start_index = np.abs(t_array[:] - t_start).argmin()
    t_stop_index  = np.abs(t_array[:] - t_stop).argmin()
    t_array = t_array[t_start_index:t_stop_index]
    fig = pl.figure(figsize=[10, 10])
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], frameon=False)
    # Plot the electrodes
    ax.scatter(Mea.elec_z, Mea.elec_y, color='b')
    # Plot the neuron
    x = np.array([cell.xstart/1000, cell.xmid/1000, cell.xend/1000])
    y = np.array([cell.ystart/1000, cell.ymid/1000, cell.yend/1000])
    z = np.array([cell.zstart/1000, cell.zmid/1000, cell.zend/1000])
    for comp in xrange(n_compartments):
        pl.plot((z[0,comp], z[2,comp]), (y[0,comp],y[2,comp]), color='k', lw=cell.diam[comp])
    signal_at_elec = signal_at_elec[:,t_start_index:t_stop_index]
    time_factor = 1./(t_stop - t_start)*Mea.electrode_separation/2
    signal_range = np.max(signal_at_elec) - np.min(signal_at_elec)
    pot_factor  = electrode_separation/signal_range
    for elec in xrange(n_elecs):
        t = t_array*time_factor + Mea.elec_z[elec]
        trace = (signal_at_elec[elec])*pot_factor + Mea.elec_y[elec]
        ax.plot(t, trace, color='r', lw = 0.5)
        #ax.plot(spike_at_t*time_factor + Mea.elec_z[elec], \
        #        np.ones(len(spike_at_t)) * Mea.elec_y[elec], 'o', color='r')
    #ax.axis([-0.1, 0.1, -0.1, 0.1])           
    #ax.axis([1.1*np.min(Mea.elec_z),1.1*np.max(Mea.elec_z) , 1.1*np.min(Mea.elec_y),1.1*np.max(Mea.elec_y) ])
    ax.axis([-.15,0.25, -0.2,0.3 ])
    ax.plot([0.075, 0.075  + electrode_separation/2], [-0.099, -0.099], color='k', lw = 3)
    ax.text(0.075, -0.105, '%g ms' % int(t_stop-t_start))
    
    ax.plot([40, 50], [-60, -60], color='k', lw = 3)
    ax.text(42, -65, '10 $\mu$m')
    
    ax.plot([0.1, 0.1], [-0.09, -0.09 + electrode_separation], color='r', lw=2)
    ax.text(0.101, -0.09 +electrode_separation/2 , '%g $\mu$V' % int(signal_range))

    pl.xlabel('z [mm]')
    pl.ylabel('y [mm]')
    #ax.set_xticks([])
    #ax.set_yticks([])

    ## ax.set_title('Location-dependent extracellular spike shapes')
    
    ## #plotting the soma trace    
    ## ax = fig.add_axes([0.75, 0.55, 0.2, 0.35])
    ## ax.plot(tvec, somav)
    ## ax.set_title('Somatic action-potential')
    ## ax.set_ylabel(r'$V_\mathrm{membrane}$ (mV)')
    ## #plotting the synaptic current
    ## ax = fig.add_axes([0.75, 0.1, 0.2, 0.35])
    ## ax.plot(tvec, synapses[0].i)
    ## ax.set_title('Synaptic current')
    ## ax.set_ylabel(r'$i_\mathrm{synapse}$ (nA)')
    ## ax.set_xlabel(r'time (ms)')
    ## #ax.axis([150,200, -1.5,1])
    #pl.savefig('morph_n_currents.png')
    pl.show()

def plot_neuron_from_side(cell, Mea, set_up_parameters):
    pl.figure(figsize=(8, 5))
    pl.subplot(121)
    n_compartments = len(cell.imem[:,0])
    #pl.plot(cell.xmid/1000, cell.ymid/1000, 'o', color='w')
    for comp in xrange(n_compartments):
        pl.plot((cell.xstart[comp]/1000, cell.xend[comp]/1000),\
                (cell.ystart[comp]/1000, cell.yend[comp]/1000),\
                color='k', lw=cell.diam[comp])
    pl.title('red=electrode plane, \n blue=tissue plane')
    #pl.axis([-0.1,1.5, -0.1, 1.5])
    pl.axis('equal')
    elec_plane = -set_up_parameters['slice_thickness']/2.
    tissue_plane = set_up_parameters['slice_thickness']/2.
    l = pl.axvline(linewidth=1, color='r', x=elec_plane)
    l = pl.axvline(linewidth=1, color='b', x=tissue_plane)
    pl.subplot(122)
    pl.axis([-0.2,0.2, -0.2, 0.6])
    pl.scatter(Mea.elec_z, Mea.elec_y, color='b')
    
    for comp in xrange(n_compartments):
        pl.plot((cell.zstart[comp]/1000, cell.zend[comp]/1000),\
                (cell.ystart[comp]/1000, cell.yend[comp]/1000),\
                color='k', lw=cell.diam[comp])
    #pl.savefig('morph_from_side.png')
    pl.show()

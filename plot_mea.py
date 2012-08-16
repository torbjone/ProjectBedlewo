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
    ax.axis('equal')
    ax.axis([-.15,0.25, -0.2,1.3 ])
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

def plot_electrodes_and_neuron(cell,  Mea):
    pl.close('all')
    n_compartments = len(cell.xstart)
    fig = pl.figure(figsize=[10, 10])
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], frameon=False)
    # Plot the electrodes
    ax.scatter(Mea.elec_z, Mea.elec_y, color='b')
    # Plot the neuron
    x = np.array([cell.xstart/1000, cell.xmid/1000, cell.xend/1000])
    y = np.array([cell.ystart/1000, cell.ymid/1000, cell.yend/1000])
    z = np.array([cell.zstart/1000, cell.zmid/1000, cell.zend/1000])
    for comp in xrange(n_compartments):
        pl.plot((z[0,comp], z[2,comp]), (y[0,comp],y[2,comp]), color='k',\
                lw=cell.diam[comp])
    ax.axis('equal')
    ax.axis([-.15,0.25, -0.2,1.3 ])
    ax.plot([40, 50], [-60, -60], color='k', lw = 3)
    ax.text(42, -65, '10 $\mu$m')
    pl.xlabel('z [mm]')
    pl.ylabel('y [mm]')
 
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
    pl.axis([-0.2,0.2, -0.2, 1.1])
    pl.scatter(Mea.elec_z, Mea.elec_y, color='b')
    
    for comp in xrange(n_compartments):
        pl.plot((cell.zstart[comp]/1000, cell.zend[comp]/1000),\
                (cell.ystart[comp]/1000, cell.yend[comp]/1000),\
                color='k', lw=cell.diam[comp])
    #pl.savefig('morph_from_side.png')
    pl.show()
    
def animate_MEA(signal, cell, Mea, plot_range):
    pl.close('all')
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    # Below are two function used in the animation. 
    def init():
        #scat = ax.scatter(x,z, c=ica[:,0], s=10*comp_size, vmin=-1, vmax=1)
        stim_point, = ax.plot([cell.zmid[stim_idx_1]], [cell.ymid[stim_idx_1]], \
                              'D', color='w', markersize=15)
        stim_point2, = ax.plot([cell.zmid[stim_idx_2]], [cell.ymid[stim_idx_2]], 'D', color='y', markersize=15)
        stim_point3, = ax.plot([cell.zmid[syn_numb]], [cell.ymid[syn_numb]], 'D', color='g', markersize=15)
        time_bar, = ax2.plot([start_t,start_t], [0,9])
        time_text.set_text('')
        return scat, time_text, stim_point
    def update_plot(i, data, scat, time_text, time_bar):
        scat.set_array(data[:,i])
        time_text.set_text(time_template%(t_array[i]))
        time_bar.set_data([t_array[i], t_array[i]], [-1000,1000])
        stim_point.set_data([cell.zmid[stim_idx_1]], [cell.ymid[stim_idx_1]])
        return scat, time_text,time_bar#, stim_point
    n_compartments = len(cell.imem[:,0])
    stim_idx_1 = 0
    stim_idx_2 = 0
    syn_numb = 0
    # Picking out the desired time range.
    start_t, stop_t = plot_range    
    start_t_ixd = np.argmin(np.abs(cell.tvec - start_t))
    stop_t_ixd  = np.argmin(np.abs(cell.tvec - stop_t))
    t_array = cell.tvec[start_t_ixd:stop_t_ixd]
    #vmem = cell.vmem[:,start_t_ixd:stop_t_ixd]
    imem = cell.imem[:,start_t_ixd:stop_t_ixd]
    signal = signal[:,start_t_ixd:stop_t_ixd]
    
    n_tsteps = len(imem[0,:])
    y = Mea.elec_y
    z = Mea.elec_z

    dt = cell.timeres_python
    # Index of the compartment recieving synaptic input or stimulation

    fig = plt.figure(figsize=[7,11])
    # Initiate cell plot with membrane voltage
    ax = fig.add_axes([0.1,0.05,0.5,0.85])
    try:
        neur_x = np.array([cell.xstart/1000, cell.xmid/1000, cell.xend/1000])
        neur_y = np.array([cell.ystart/1000, cell.ymid/1000, cell.yend/1000])
        neur_z = np.array([cell.zstart/1000, cell.zmid/1000, cell.zend/1000])
        for comp in xrange(n_compartments):
            pl.plot((neur_z[0,comp], neur_z[2,comp]), \
                    (neur_y[0,comp],neur_y[2,comp]), color='k', lw=cell.diam[comp])
    except:
        for comp in xrange(n_compartments):
            pl.plot(cell.zmid, cell.ymid, 'o')
    #scat = ax.scatter(x,z, c=ica[:,0], s=10*comp_size, vmin=-0.01, vmax=0.01)
    scat = ax.scatter(z,y, c=signal[:,0], s=500, vmin = np.min(signal), vmax = np.max(signal), marker = 'D')
    stim_point, = ax.plot([cell.zmid[stim_idx_1]], [cell.ymid[stim_idx_1]], 'D', color='w')
    ax.axis('equal')
    ax.axis([-0.05, .1, .500,1.0])
    pl.colorbar(scat)

    # Plot of soma membrane current
    ax2 = fig.add_axes([0.7,0.8,0.25,0.15])
    ax2.set_title('imem at soma')
    ax2.axis([start_t, stop_t, np.min(imem[0,:]), np.max(imem[0,:])])
    ax2.plot(t_array, imem[0,:])
    time_bar, = ax2.plot([start_t,start_t], [0,9])
    time_template = 'time = %.3fms'
    time_text = ax.text(0, max(z)*1.1, '')

    # Stimulation point-1 current
    ax3 = fig.add_axes([0.7,0.6,0.25,0.15])
    ax3.set_title('imem at inj. point 1 (wh)')
    ax3.axis([start_t, stop_t, np.min(imem[stim_idx_1,:]), np.max(imem[stim_idx_1,:])])
    ax3.plot(t_array, imem[stim_idx_1,:])

    # Stimulation point-2 current
    ax4 = fig.add_axes([0.7,0.4,0.25,0.15])
    ax4.set_title('imem at inj. point 2 (ye)')
    ax4.axis([start_t, stop_t, np.min(imem[stim_idx_2,:]), np.max(imem[stim_idx_2,:])])
    ax4.plot(t_array, imem[stim_idx_2,:])

    # Chosen synaps membrane current
    ax4 = fig.add_axes([0.7,0.2,0.25,0.15])
    ax4.set_title('mem at point %d (gr)' %syn_numb)
    ax.plot([cell.xmid[syn_numb]], [cell.zmid[syn_numb]], 'D', color='g')
    ax4.axis([start_t, stop_t, np.min(imem[syn_numb,:]), np.max(imem[syn_numb,:])])
    ax4.plot(t_array, imem[syn_numb,:])

    ani = animation.FuncAnimation(fig, update_plot, frames=xrange(n_tsteps),
                                  fargs=(signal, scat, time_text, time_bar), \
                                  blit=True, interval=50, init_func=init)
    #ani.save('simple_plot.mp4')
    #ani.ffmpeg_cmd('simple_plot.mp4', fps=5, codec='mpeg4',  frame_prefix='_tmp')    
    pl.show()


import matplotlib.pylab as pl
import numpy as np
import LFPy
import os
#from mayavi import mlab
import neuron
pl.rcParams.update({'font.size' : 12,
    'figure.facecolor' : '1',
    'wspace' : 0.5, 'hspace' : 0.5})

#seed for random generation
pl.seed(0)

def plotstuff(cell):
    fig = pl.figure(figsize=[12, 8])
    
    ax = fig.add_axes([0.1, 0.7, 0.5, 0.2])
    ax.plot(cell.tvec,cell.imem[0,:])
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('Soma trans. current')
    
    ax = fig.add_axes([0.1, 0.4, 0.5, 0.2])
    for i in xrange(len(cell.synapses)):
        ax.plot(cell.tvec,cell.synapses[i].i,color=cell.synapses[i].color)
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('Syn. i [nA]')

    stim_array = np.sum(cell.imem, axis=0)
    ax = fig.add_axes([0.1, 0.1, 0.5, 0.2])
    ax.plot(cell.tvec,stim_array)
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('Syn. i [nA]')
    ax = fig.add_axes([0.65, 0.1, 0.25, 0.8], frameon=False)
    for i in xrange(cell.xend.size):
        ax.plot([cell.xstart[i],cell.xend[i]],[cell.zstart[i],cell.zend[i]],color='k')
    for i in xrange(len(cell.synapses)):
        ax.plot([cell.synapses[i].x],[cell.synapses[i].z],\
            color=cell.synapses[i].color,marker=cell.synapses[i].marker)
    stim_idx_1 = clamp_1['idx']
    stim_idx_2= clamp_2['idx']
    ax.plot([cell.xmid[stim_idx_1]], [cell.zmid[stim_idx_1]], 'D',color='y')
    ax.plot([cell.xmid[stim_idx_2]], [cell.zmid[stim_idx_2]], 'D',color='g')
    pl.axis('equal')
    pl.axis(pl.array(pl.axis())*0.8)
    ax.set_xticks([])
    ax.set_yticks([])
    pl.savefig('example_fig.png')

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
    'tstopms' : 1000,           #stop time of simulation
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
    'n' : 10,
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
    'idx' : 221,
    'record_current' : True,
    'amp' : 3., #[nA]
    'dur' : 500,
    'delay' :100,
    #'freq' : 10,
    #'phase' : 0,
    #'pkamp' : 300e-3,
    'pptype' : 'IClamp',
}

clamp_2 = {
    'idx' : 375,
    'record_current' : True,
    'amp' : 3., #[nA]
    'dur' : 500,
    'delay' :100,
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
        idx = cell.get_rand_idx_area_norm(section=section, nidx=n)

        #Insert synapses in an iterative fashion
        for i in idx:
            synparams.update({'idx' : int(i)})

            # Some input spike train using the function call
            spiketimes = spTimesFun(args[0], args[1], args[2], args[3])

            # Create synapse(s) and setting times using the Synapse class in LFPy
            s = LFPy.Synapse(cell,**synparams)
            s.set_spike_times(spiketimes)

    cell = LFPy.Cell(**cellParameters)
    cell.set_rotation(x = pl.pi/2)
    #cell.set_rotation(z = pl.pi/2)
    if do_simulation:
        np.save(output_folder + 'x_start.npy', cell.xstart)
        np.save(output_folder + 'y_start.npy', cell.ystart)
        np.save(output_folder + 'z_start.npy', cell.zstart)
        np.save(output_folder + 'x_end.npy', cell.xend)
        np.save(output_folder + 'y_end.npy', cell.yend)
        np.save(output_folder + 'z_end.npy', cell.zend)
        np.save(output_folder + 'diam.npy', cell.diam)
        currentClamp_1 = LFPy.StimIntElectrode(cell, **clamp_1)
        currentClamp_2 = LFPy.StimIntElectrode(cell, **clamp_2)
        #insert_synapses(synapseParameters_AMPA, **insert_synapses_AMPA_args)
        #insert_synapses(synapseParameters_NMDA, **insert_synapses_NMDA_args)
        #insert_synapses(synapseParameters_GABA_A, **insert_synapses_GABA_A_args)
        cell.simulate(**simulationParameters)
        np.save(output_folder + 'imem.npy', cell.imem)
        np.save(output_folder + 'vmem.npy', cell.vmem)
        np.save(output_folder + 'tvec.npy', cell.tvec)
        plotstuff(cell)
    else:
        cell.vmem = np.load(output_folder + 'vmem.npy')
        cell.imem = np.load(output_folder + 'imem.npy')
        cell.tvec = np.load(output_folder + 'tvec.npy')
    return cell

def plot_cell_3D(cell):
    savefolder = 'savefiles'
    ## #######MAYAVISTUFF###############################################
    from visio import VisioAnim as Visio
    import scitools

    allsecnames = []
    allsec = []
    for sec in cell.allseclist:
        allsecnames.append(sec.name())
        for seg in sec:
            allsec.append(sec.name())
    cell.allsecnames = allsecnames
    cell.allsec = allsec

    visCurrent2 = Visio()
    try:
        os.mkdir('mkdir %s' % savefolder + '/anim_imem2')
        os.system('rm -r %s' % savefolder + '/anim_imem2/*')
    except:
        pass
    visCurrent2.draw_stick_model(cell, values = 1000*(cell.imem.T / cell.area).T)
    for v in visCurrent2.voltages:
        visCurrent2.draw_surface(0, v[:, 0], 'v',
                                colormap='jet',
                                vmin=-1, vmax=1)

        visCurrent2.mayavi.visualization.scene.disable_render = True    
    #change background color, orientation, create colorbar
    visCurrent2.mayavi.visualization.scene.background = (0,0,0)
    visCurrent2.mayavi.visualization.scene.mlab.view(0,0, 2500)
    visCurrent2.mayavi.visualization.scene.mlab.colorbar(orientation='vertical',
                        title='i(pA/mum2)', object=visCurrent2.surf, nb_labels=11)
    
    
    #speed up?
    visCurrent2.mayavi.visualization.scene.disable_render = True
        
    #update membrane voltages for every timestep
    for j in xrange(60,4*60):
        i = 4*j #plot every n timestep
        print '.',
        #draw the morphologies
  
        visCurrent2.draw_surface(0,
                                visCurrent2.voltages[0][:, i], 'v',
                                colormap='jet',
                                    vmin=-1, vmax=1)
        
        #set the colorbar
        visCurrent2.mayavi.visualization.scene.mlab.colorbar(
                orientation='vertical', title='i(pA/mum2)',
                object=visCurrent2.surf, nb_labels=11)
        
        #set the view
        visCurrent2.mayavi.visualization.scene.mlab.view(0,90, 2500)
        #display the time
        visCurrent2.mayavi.visualization.scene.mlab.title(
            't = %.f ms' % np.floor(i * cell.timeres_python))
        visCurrent2.mayavi.visualization.scene.mlab.savefig(
            savefolder + '/anim_imem2/cell_tstep_%.5i.png' % j)
    visCurrent2.mayavi.visualization.scene.disable_render = True
        
    #compile into movie
    #scitools.easyviz.movie(savefolder + '/anim_imem2/cell_tstep_*.png',
    #                       encoder='ffmpeg', output_file='movie_imem2.avi',
    #                       vcodec='mjpeg',qscale=5)

    os.system('mencoder "mf://%s/anim_imem2/*.png" -mf type=png:fps=10 -ovc lavc -o output.avi'% savefolder)

def simple_plot_2D(cell, plot_range):
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    def init():
        ax.plot([cell.xmid[stim_idx_1]], [cell.zmid[stim_idx_2]], 'D', color='y')
        scat = ax.scatter(x,z, c=imem[:,0], s=10*comp_size)
        stim_point, = ax.plot([cell.xmid[stim_idx_1]], [cell.zmid[stim_idx_1]], 'D', color='w')
        time_bar, = ax2.plot([start_t,start_t], [0,9])
        time_text.set_text('')
        return scat, time_text, stim_point
    def update_plot(i, data, scat, time_text, time_bar, stim_point):
        scat.set_array(data[:,i])
        time_text.set_text(time_template%(t_array[i]))
        time_bar.set_data([t_array[i], t_array[i]], [-1000,1000])
        stim_point.set_data([cell.xmid[stim_idx_1]], [cell.zmid[stim_idx_1]])
        return scat, time_text,time_bar, stim_point
    start_t, stop_t = plot_range    
    start_t_ixd = np.argmin(np.abs(cell.tvec - start_t))
    stop_t_ixd  = np.argmin(np.abs(cell.tvec - stop_t))
    t_array = cell.tvec[start_t_ixd:stop_t_ixd]
    imem = cell.vmem[:,start_t_ixd:stop_t_ixd]
    n_tsteps = len(imem[0,:])
    x = cell.xmid
    y = cell.ymid
    z = cell.zmid
    dt = cell.timeres_python
    comp_size = cell.diam
    stim_idx_1 = clamp_1['idx']
    stim_idx_2 = clamp_2['idx']
    
    pl.close('all')
    fig = plt.figure(figsize=[7,11])
    ax = fig.add_axes([0.1,0.1,0.5,0.9])
    scat = ax.scatter(x,z, c=imem[:,0], s=10*comp_size)
    stim_point, = ax.plot([cell.xmid[stim_idx_1]], [cell.zmid[stim_idx_1]], 'D', color='w')
    #ax.axis('equal')
    ax.axis([-100, 300, -200,1100])
    pl.colorbar(scat)

    # Soma current
    ax2 = fig.add_axes([0.7,0.8,0.25,0.15])
    ax2.set_title('vmem at soma')
    ax2.axis([start_t, stop_t, np.min(imem[0,:]), np.max(imem[0,:])])
    stim_point, = ax.plot([cell.xmid[stim_idx_1]], [cell.zmid[stim_idx_1]], 'D', color='w')
    ax2.plot(t_array, imem[0,:])
    time_bar, = ax2.plot([start_t,start_t], [0,9])
    time_template = 'time = %.3fms'
    time_text = ax.text(0, max(z)*1.1, '')

    # Stimulation point-1 current
    ax3 = fig.add_axes([0.7,0.6,0.25,0.15])
    ax3.set_title('vmem at inj. point 1')
    ax3.axis([start_t, stop_t, np.min(imem[stim_idx_1,:]), np.max(imem[stim_idx_1,:])])
    ax3.plot(t_array, imem[stim_idx_1,:])

    # Stimulation point-2 current
    ax4 = fig.add_axes([0.7,0.4,0.25,0.15])
    ax4.set_title('vmem at inj. point 2')
    ax4.axis([start_t, stop_t, np.min(imem[stim_idx_2,:]), np.max(imem[stim_idx_2,:])])
    ax4.plot(t_array, imem[stim_idx_2,:])

    ani = animation.FuncAnimation(fig, update_plot, frames=xrange(n_tsteps),
                                  fargs=(imem, scat, time_text, time_bar, stim_point), blit=True, interval=.1, init_func=init)
    #ani.save('simple_plot.mp4')
    #ani.ffmpeg_cmd('simple_plot.mp4', fps=5, codec='mpeg4',  frame_prefix='_tmp')    
    pl.show()

def plot_cell_compartments(cell):
    for comp_idx in xrange(len(cell.xmid)):
        pl.plot(cell.xmid[comp_idx], cell.zmid[comp_idx], marker='$%i$'%comp_idx, color='b', markersize=10)
    pl.show()

if __name__ == '__main__':
    output_folder = 'larkum_sim/'
    do_simulation = True
    plot_range = [95,200]
    try:
        os.mkdir(output_folder)
    except(OSError):
        if do_simulation:
            print "Result folder already exists. Overwriting..."
        else:
            print "Loading simulation files..."
    cell = get_cell(output_folder, do_simulation)
    #plot_cell_compartments(cell)
    simple_plot_2D(cell, plot_range)

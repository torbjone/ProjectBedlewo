
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
pl.seed(999)

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
    
   
    ax = fig.add_axes([0.65, 0.1, 0.25, 0.8], frameon=False)
    for i in xrange(cell.xend.size):
        ax.plot([cell.xstart[i],cell.xend[i]],[cell.zstart[i],cell.zend[i]],color='k')
    for i in xrange(len(cell.synapses)):
        ax.plot([cell.synapses[i].x],[cell.synapses[i].z],\
            color=cell.synapses[i].color,marker=cell.synapses[i].marker)
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
    'timeres_NEURON' : 2**-4,   # dt of LFP and NEURON simulation.
    'timeres_python' : 2**-4,
    'tstartms' : -10,          #start time, recorders start at t=0
    'tstopms' : 100,           #stop time of simulation
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


clampparams = {
    'idx' : 0,
    'record_current' : True,
    'amp' : 400e-2, #[mA]
    'dur' : 50,
    'delay' :4,
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

def get_cell(do_simulation = True):
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
    cell.set_rotation(z = pl.pi/2)
    if do_simulation:
        np.save('x_start.npy', cell.xstart)
        np.save('y_start.npy', cell.ystart)
        np.save('z_start.npy', cell.zstart)
        np.save('x_end.npy', cell.xend)
        np.save('y_end.npy', cell.yend)
        np.save('z_end.npy', cell.zend)
        np.save('diam.npy', cell.diam)
        currentClamp = LFPy.StimIntElectrode(cell, **clampparams)
        insert_synapses(synapseParameters_AMPA, **insert_synapses_AMPA_args)
        insert_synapses(synapseParameters_NMDA, **insert_synapses_NMDA_args)
        insert_synapses(synapseParameters_GABA_A, **insert_synapses_GABA_A_args)
        cell.simulate(**simulationParameters)
        np.save('imem.npy', cell.imem)
        np.save('tvec.npy', cell.tvec)
        plotstuff(cell)
    else:
        cell.imem = np.load('imem.npy')
        cell.tvec = np.load('tvec.npy')
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

def simple_plot_2D(cell, start_t, stop_t):
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    def init():
        scat = ax.scatter(y,z, c=imem[:,0], s=10*comp_size)
        ax2.plot(t_array, imem[0,:])
        time_bar = ax2.plot([start_t,start_t], [0,9])
        time_text.set_text('')
        return scat, time_text
    def update_plot(i, data, scat, time_bar):
        scat.set_array(data[:,i])
        time_text.set_text(time_template%(t_array[i]))
        time_bar.set_data([t_array[i], t_array[i]], [0,9])
        return scat, time_text,time_bar
   
    start_t_ixd = np.argmin(np.abs(cell.tvec - start_t))
    stop_t_ixd  = np.argmin(np.abs(cell.tvec - stop_t))
    t_array = cell.tvec[start_t_ixd:stop_t_ixd]
    imem = cell.imem[:,start_t_ixd:stop_t_ixd]
    n_tsteps = len(imem[0,:])
    x = cell.xmid
    y = cell.ymid
    z = cell.zmid
    dt = cell.timeres_python
    comp_size = cell.diam
    fig = plt.figure(figsize=[7,11])
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax2.axis([start_t, stop_t, 0,9])
    #ax2.plot(t_array, imem[0,:])
    time_bar, = ax2.plot([start_t,start_t], [0,9])
    time_template = 'time = %.3fms'
    time_text = ax.text(0, max(z)*1.1, '')
    scat = ax.scatter(y,z, c=imem[:,0], s=10*comp_size)
    ax.axis('equal')
    pl.colorbar(scat)
    ani = animation.FuncAnimation(fig, update_plot, frames=xrange(n_tsteps),
                                  fargs=(imem, scat, time_bar),blit=True, interval=.01, init_func=init)
    #ani.save('simple_plot.mp4')
    #ani.ffmpeg_cmd('simple_plot.mp4', fps=5, codec='mpeg4',  frame_prefix='_tmp')    
    pl.show()

if __name__ == '__main__':
    cell = get_cell(do_simulation = False)
    simple_plot_2D(cell, 35,42)

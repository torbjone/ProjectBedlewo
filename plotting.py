
import matplotlib.pylab as pl
import numpy as np
import LFPy
import os
import sys
import neuron
pl.rcParams.update({'font.size' : 12,
    'figure.facecolor' : '1',
    'wspace' : 0.5, 'hspace' : 0.5})

def plot_cell_compartments(cell):
    for comp_idx in xrange(len(cell.xmid)):
        pl.plot(cell.zmid[comp_idx], cell.ymid[comp_idx],\
                marker='$%i$'%comp_idx, color='b', markersize=10)
    pl.show()
    pl.close('all')

def plotstuff(cell, clamp_1, clamp_2):
    fig = pl.figure(figsize=[12, 8])
    
    ax = fig.add_axes([0.1, 0.7, 0.5, 0.2])
    ax.plot(cell.tvec,cell.vmem[0,:])
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
    ax = fig.add_axes([0.65, 0.05, 0.25, 1], frameon=False)
    for i in xrange(cell.xend.size):
        ax.plot([cell.zstart[i],cell.zend[i]],[cell.ystart[i],cell.yend[i]],color='k')
    for i in xrange(len(cell.synapses)):
        ax.plot([cell.synapses[i].z],[cell.synapses[i].y],\
            color=cell.synapses[i].color,marker=cell.synapses[i].marker)
    stim_idx_1 = clamp_1['idx']
    stim_idx_2 = clamp_2['idx']
    ax.plot([cell.zmid[stim_idx_1]], [cell.ymid[stim_idx_1]], 'D',color='y')
    ax.plot([cell.zmid[stim_idx_2]], [cell.ymid[stim_idx_2]], 'D',color='g')
    pl.axis('equal')
    pl.axis(pl.array(pl.axis())*1.2)
    ax.set_xticks([])
    ax.set_yticks([])
    pl.savefig('example_fig.png')
    pl.close('all')

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

def simple_plot_2D(cell, plot_range, clamp_1, clamp_2):
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

    # Picking out the desired time range.
    start_t, stop_t = plot_range    
    start_t_ixd = np.argmin(np.abs(cell.tvec - start_t))
    stop_t_ixd  = np.argmin(np.abs(cell.tvec - stop_t))
    t_array = cell.tvec[start_t_ixd:stop_t_ixd]
    
    vmem = cell.vmem[:,start_t_ixd:stop_t_ixd]
    #ica = cell.rec_variables['ica'][:,start_t_ixd:stop_t_ixd]
    #cai = cell.rec_variables['cai'][:,start_t_ixd:stop_t_ixd]
    imem = cell.imem[:,start_t_ixd:stop_t_ixd]
    signal = vmem
    n_tsteps = len(vmem[0,:])
    x = cell.xmid
    y = cell.ymid
    z = cell.zmid
    dt = cell.timeres_python
    comp_size = cell.diam
    # Index of the compartment recieving synaptic input or stimulation
    stim_idx_1 = clamp_1['idx']
    stim_idx_2 = clamp_2['idx']
    
    fig = plt.figure(figsize=[7,11])
    # Initiate cell plot with membrane voltage
    ax = fig.add_axes([0.1,0.1,0.5,0.9])
    #scat = ax.scatter(x,z, c=ica[:,0], s=10*comp_size, vmin=-0.01, vmax=0.01)
    scat = ax.scatter(z,y, c=signal[:,0], s=10*comp_size, \
                      vmin=np.min(signal), vmax=np.max(signal))
    stim_point, = ax.plot([cell.zmid[stim_idx_1]], [cell.ymid[stim_idx_1]],\
                          'D', color='w')
    #ax.axis('equal')
    ax.axis([-200, 300, -300,1200])
    pl.colorbar(scat)

    # Plot of soma membrane current
    ax2 = fig.add_axes([0.7,0.8,0.25,0.15])
    ax2.set_title('soma')
    ax2.axis([start_t, stop_t, np.min(signal[0,:]), np.max(signal[0,:])])
    ax2.plot(t_array, signal[0,:])
    time_bar, = ax2.plot([start_t,start_t], [0,9])
    time_template = 'time = %.3fms'
    time_text = ax.text(0, max(z)*1.1, '')

    # Stimulation point-1 current
    ax3 = fig.add_axes([0.7,0.6,0.25,0.15])
    ax3.set_title('inj. point 1 (wh)')
    ax3.axis([start_t, stop_t, np.min(signal[stim_idx_1,:]), \
              np.max(signal[stim_idx_1,:])])
    ax3.plot(t_array, signal[stim_idx_1,:])

    # Stimulation point-2 current
    ax4 = fig.add_axes([0.7,0.4,0.25,0.15])
    ax4.set_title('inj. point 2 (ye)')
    ax4.axis([start_t, stop_t, np.min(signal[stim_idx_2,:]), \
              np.max(signal[stim_idx_2,:])])
    ax4.plot(t_array, signal[stim_idx_2,:])

    # Chosen synaps membrane current
    syn_numb = 179
    ax4 = fig.add_axes([0.7,0.2,0.25,0.15])
    ax4.set_title(' point %d (gr)' %syn_numb)
    ax.plot([cell.xmid[syn_numb]], [cell.zmid[syn_numb]], 'D', color='g')
    ax4.axis([start_t, stop_t, np.min(signal[syn_numb,:]), \
              np.max(signal[syn_numb,:])])
    ax4.plot(t_array, signal[syn_numb,:])

    ani = animation.FuncAnimation(fig, update_plot, frames=xrange(n_tsteps),
                                  fargs=(signal, scat, time_text, time_bar), blit=True, interval=.1, init_func=init)
    #ani.save('simple_plot.mp4')
    #ani.ffmpeg_cmd('simple_plot.mp4', fps=5, codec='mpeg4',  frame_prefix='_tmp')    
    pl.show()

#!/usr/bin/env python
import os
os.environ['ETS_TOOLKIT'] = 'qt4'
import logging
logger = logging.getLogger(__name__)

from pyface.qt import QtGui, QtCore
from PyQt4 import uic

from traits.api import HasTraits, Instance, on_trait_change, \
    Int, Dict
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
        SceneEditor

import numpy as np
from mayavi import mlab
#mlab.options.backend = 'envisage'

from neuron import h



################################################################################
#The actual visualization
class Visualization(HasTraits):
    
    scene = Instance(MlabSceneModel, ())


    @on_trait_change('scene.activated')
    def update_plot(self):
        # This function is called when the view is opened. We don't
        # populate the scene when the view is not yet open, as some 
        # VTK features require a GLContext.

        # We can do normal mlab calls on the embedded scene.
#        self.scene.mlab.test_points3d()
        pass

    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=1000, width=1000,
                     show_label=False),
                resizable=True, # We need this to resize with the parent widget
                )


################################################################################
# The QWidget containing the visualization, this is pure PyQt4 code.
class MayaviQWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()
        
        # If you want to debug, beware that you need to remove the Qt
        # input hook.
        #QtCore.pyqtRemoveInputHook()
        #import pdb ; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self, 
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)



class Visio(object):
    
    def __init__(self):       
        
        ## Needed when user pick the cylinder from visio and 
        ## we need to get the section
        self.cyl2sec = {}
        
        ## Needed to update the value of a cyl bound to a section
        self.sec2cyl = {}
        
        self.seg2id = {}
        self.sec2coords = {}
        #self.connections = []
        self.n3dpoints_per_sec = {}
           
        self.selected_cyl = None # Used for storing the cyl when picked
        #self.sec_info_label = sec_info_label # Info for the selected sec
        #self.manager = manager
        #self.manager = {'v'}
        
        container = QtGui.QWidget()
        container.setWindowTitle("Hell Yeah!")
        
        self.mayavi = MayaviQWidget(container)
        layout = QtGui.QVBoxLayout(container)
        layout.addWidget(self.mayavi)
                
        # Tell visual to use this as the viewer.
        #visual.set_viewer(self.mayavi.visualization.scene)
        
        # binding to hide event.
        container.connect(container, QtCore.SIGNAL('closeEvent()'), 
                               self.closeEvent)
        
        container.show()
        
        self.container = container
        
        ## Connecting the picker.
        #figure = self.mayavi.visualization.scene.mlab.gcf()
        #self.outline = None
        #self.picker = figure.on_mouse_pick(self.picker_callback, type='cell')
        #
        ## ScalarBar and time_point
        #self.colorbar = None
        #self.timelabel = None
    
    def get_idx_section(self, allsec, section='soma[0]'):
        idxvec = np.zeros(len(allsec))
        secnamelist = []
        i = 0
        for sec in allsec:
            if sec == section:
                idxvec[i] += 1
            i += 1
        [idx] = np.where(idxvec)
        return idx
    
    
    def draw_stick_model(self, cell, color=(1,0,0)):
        '''draw the model as set of line sources'''
        # Disable the render. Faster drawing.
        self.mayavi.visualization.scene.disable_render = True
        

        x,y,z,d = [], [], [], []
        voltage = []
        connections = []
        
        for sec in cell.allsecnames:
            #section = sec.name()
            idx = self.get_idx_section(cell.allsec, section=sec)
            for i, xi in enumerate(idx):
                x.append(cell.xstart[idx[i]])
                y.append(cell.ystart[idx[i]])
                z.append(cell.zstart[idx[i]])
                d.append(cell.diam[idx[i]])
                indx_geom_seg = len(x) - 1
                if len(x) > 1 and i > 0:
                    connections.append([indx_geom_seg, indx_geom_seg-1])
            x.append(cell.xend[idx[i]])
            y.append(cell.yend[idx[i]])
            z.append(cell.zend[idx[i]])
            d.append(cell.diam[idx[i]])
            indx_geom_seg = len(x) - 1
            connections.append([indx_geom_seg, indx_geom_seg-1])
            
        self.edges  = connections
        self.x = x
        self.y = y
        self.z = z
        
        # Mayavi pipeline        
        d = np.array(d) # Transforming for easy division

        self.draw_mayavi(x, y, z, d, self.edges, color)        
        
        

    def draw_model(self):
        """Draw the model.
        Params:
        controls - the main gui obj."""
        
        # Draw the new one
        h.define_shape()
        num_sections = 0

        # Disable the render. Faster drawing.
        self.mayavi.visualization.scene.disable_render = True
        

        x,y,z,d = [], [], [], []
        voltage = []
        connections = []
        for sec in h.allsec():
            x_sec, y_sec, z_sec, d_sec = self.retrieve_coordinate(sec)
            self.sec2coords[sec.name()] = [x_sec, y_sec, z_sec]
            # Store the section. later.
            radius = sec.diam/2.
            sec_coords_bound = ((x_sec.min(), x_sec.max()), 
                                (y_sec.min() - radius, 
                                 y_sec.max() + radius), 
                                (z_sec.min() - radius, 
                                 z_sec.max() + radius))
            self.cyl2sec[sec_coords_bound] = sec 
            self.sec2cyl[sec] = sec_coords_bound
            
            
            for i,xi in enumerate(x_sec):
                x.append(x_sec[i])
                y.append(y_sec[i])
                z.append(z_sec[i])
                d.append(d_sec[i])
                indx_geom_seg = len(x) -1
                
                if len(x) > 1 and i > 0:
                    connections.append([indx_geom_seg, indx_geom_seg-1])
            
        self.edges  = connections
        self.x = x
        self.y = y
        self.z = z
        
        # Mayavi pipeline        
        d = np.array(d) # Transforming for easy division        
        
        self.draw_mayavi(x, y, z, d, self.edges)

    def draw_mayavi(self, x, y, z, d, edges, color=(1, 0, 0)):
        "Draw the surface the first time"
        
        # rendering disabled
        self.mayavi.visualization.scene.disable_render = True
        
        points = mlab.pipeline.scalar_scatter(x, y, z, d/2.0)
        dataset = points.mlab_source.dataset
        dataset.point_data.get_array(0).name = 'diameter'
        dataset.lines = np.vstack(edges)
        dataset.point_data.update()
        self.dataset = dataset

        # The tube
        src = mlab.pipeline.set_active_attribute(points, point_scalars='diameter')
        stripper = mlab.pipeline.stripper(src)
        tube = mlab.pipeline.tube(stripper, tube_sides = 6, tube_radius = 1)
        tube.filter.capping = True
#        tube.filter.use_default_normal = False
        tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
        self.tube = tube
        

        # Setting the voltage
        # Making room for the voltage
        v = []
        for sec in h.allsec():
            sec.push()
            v.extend(np.repeat(0.0, h.n3d()))
            h.pop_section()
        
        v = np.array(v)
        
        v = np.ones(len(x)) * np.random.rand()
        
        self.draw_surface(v, 'v', color)
        
        # ReEnable the rendering
        self.mayavi.visualization.scene.disable_render = True

        
    def draw_surface(self, scalar, scalar_name, color):
        
        self.tube.children= [] # Removing the old ones 
        
        #scalar = self.get_var_data('v', 0)
        scalar = np.array(scalar)
        array_id = self.dataset.point_data.add_array(scalar)
        self.dataset.point_data.get_array(array_id).name = scalar_name
        self.dataset.point_data.update()
        
        src2 = mlab.pipeline.set_active_attribute(self.tube, 
                                                  point_scalars=scalar_name)
        #self.surf = mlab.pipeline.surface(src2, colormap='blue-red')
        self.surf = mlab.pipeline.surface(src2, color=color)
    

    def closeEvent(self):
        """Just hide the window to not loose the mayavi hook"""
        self.container.hide()

    def retrieve_coordinate(self, sec):
        """Retrieve the coordinates of the section avoiding duplicates"""
        
        sec.push()
        x, y, z, d = [],[],[],[]

        tot_points = 0
        connect_next = False
        for i in range(int(h.n3d())):
            present = False
            x_i = h.x3d(i)
            y_i = h.y3d(i)
            z_i = h.z3d(i)
            d_i = h.diam3d(i)
            # Avoiding duplicates in the sec
            if x_i in x:
                ind = len(x) - 1 - x[::-1].index(x_i) # Getting the index of last value
                if y_i == y[ind]:
                    if z_i == z[ind]:
                        present = True
                    
            if not present:
                k =(x_i, y_i, z_i)
                x.append(x_i)
                y.append(y_i)
                z.append(z_i)
                d.append(d_i)                
        h.pop_section()
        #adding num 3d points per section
        self.n3dpoints_per_sec[sec.name()] = len(d)
        return (np.array(x),np.array(y),np.array(z),np.array(d))


    def get_var_data(self, var, time_point=0):
        """Retrieve the value of the `var` for the `time_point`.
        Prameters:
        var - variable to retrieve
        time_point - point in the simulation"""
        
        var_scalar = []
        for sec in h.allsec():
            var_value = 0
            #if self.manager.refs.has_key('VecRef'):
            #    for vecRef in self.manager.refs['VecRef']:
            #        if vecRef.sec.name() == sec.name():
            #            if vecRef.vecs.has_key(var):
            #                vec = vecRef.vecs[var]
            #                try:
            #                    var_value = vec[time_point]
            #                except IndexError:
            #                    pass # vector exist, but not initialized.
            sec_scalar = self.build_sec_scalar(sec, var_value)
            var_scalar.extend(sec_scalar)
                
                    
                        
        
        if len(var_scalar) == 0:
            logger.debug( "Var scalar 0 length. Var: %s point_time: %s" %(var, 
                                                                  time_point))
        return np.array(var_scalar)


    def build_sec_scalar(self, sec, var_value):
        
        sec.push()
        npoints = self.n3dpoints_per_sec[sec.name()]
        sec_scalar = np.repeat(var_value, npoints)
        h.pop_section()
        return sec_scalar


class VisioAnim(Visio):
    '''subclassing Visio, as animatings are broken and will not break old usage'''
    
    def __init__(self):
        Visio.__init__(self)
        
    def draw_stick_model(self, cell, values='vmem'):
        '''draw the model as set of line sources'''
        # Disable the render. Faster drawing.
        self.mayavi.visualization.scene.disable_render = True
        x,y,z,d = [], [], [], []
        v = []
        connections = []
        
        for sec in cell.allsecnames:
            #section = sec.name()
            idx = self.get_idx_section(cell.allsec, section=sec)
            for i, xi in enumerate(idx):
                x.append(cell.xstart[idx[i]])
                y.append(cell.ystart[idx[i]])
                z.append(cell.zstart[idx[i]])
                d.append(cell.diam[idx[i]])
                v.append(values[idx[i], :])
                indx_geom_seg = len(x) - 1
                if len(x) > 1 and i > 0:
                    connections.append([indx_geom_seg, indx_geom_seg-1])
            x.append(cell.xend[idx[i]])
            y.append(cell.yend[idx[i]])
            z.append(cell.zend[idx[i]])
            d.append(cell.diam[idx[i]])
            v.append(values[idx[i], :])
            indx_geom_seg = len(x) - 1
            connections.append([indx_geom_seg, indx_geom_seg-1])

        # Mayavi pipeline        
        d = np.array(d) # Transforming for easy division
        
        #self.draw_mayavi(x, y, z, d, v, self.edges)      
        self.draw_mayavi(x, y, z, d, v, connections)      

        # ReEnable the rendering
        self.mayavi.visualization.scene.disable_render = True


    def draw_mayavi(self, x, y, z, d, v, edges):
        "Draw the surface the first time"
        
        # rendering disabled
        self.mayavi.visualization.scene.disable_render = True
        
        points = mlab.pipeline.scalar_scatter(x, y, z, d)
        dataset = points.mlab_source.dataset
        dataset.point_data.get_array(0).name = 'diameter'
        dataset.lines = np.vstack(edges)
        dataset.point_data.update()
        if not hasattr(self, 'datasets'):            
            self.datasets = []
        self.datasets.append(dataset)

        # The tube
        src = mlab.pipeline.set_active_attribute(points,
                                                 point_scalars='diameter')
        stripper = mlab.pipeline.stripper(src)
        tube = mlab.pipeline.tube(stripper, tube_sides = 6, tube_radius = 1)
        tube.filter.capping = True
#        tube.filter.use_default_normal = False
        tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
        if not hasattr(self, 'tubes'):
            self.tubes = []
        self.tubes.append(tube)
        

        # Setting the voltage
        # Making room for the voltage
        if not hasattr(self, 'voltages'):
            self.voltages = []
        self.voltages.append(np.array(v))
    
        ##draw the surface
        #self.draw_surface(self.v[:, 0], 'v')
        #
        ## ReEnable the rendering
        #self.mayavi.visualization.scene.disable_render = False

    def draw_surface(self, i, scalar, scalar_name, colormap='spectral', vmin=-80., vmax=20.):
        
        tube = self.tubes[i]
        tube.children= [] # Removing the old ones 
        
        #scalar = self.get_var_data('v', 0)
        scalar = np.array(scalar)
        #for dataset in self.datasets:
        dataset = self.datasets[i]
        array_id = dataset.point_data.add_array(scalar)
        dataset.point_data.get_array(array_id).name = scalar_name
        dataset.point_data.update()
        
        src2 = mlab.pipeline.set_active_attribute(tube, 
                                                  point_scalars=scalar_name)
        self.surf = mlab.pipeline.surface(src2, opacity=1.,
                                          colormap=colormap,
                                          vmin=vmin, vmax=vmax)
    

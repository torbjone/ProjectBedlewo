"""SpikeSorting module. Contains the spike sorting classes with parameters."""
import numpy as np

class HD_MEA_CMOS128:
    """HD_MEA_CMOS128 class. High-density Multielectrode Array with highest posible density. Honeycomb patters"""
    def __init__(self, elec_x):
        self.pattern = 'honeycomp'
        self.electrode_pitch = 0.017
        self.electrode_separation = self.electrode_pitch*2 #0.2/100
        self.n_elec_rows = 1
        self.n_elec_cols = 1
        self.n_elecs = self.n_elec_cols * self.n_elec_rows
        self.cylinder_height = 3.
        self.sigma_E = .3
        self.sigma_S = 3.
        self.elec_x = elec_x
        self.elec_z = np.zeros(self.n_elecs)
        self.elec_y = np.zeros(self.n_elecs)
        i = 0
        for col in xrange(self.n_elec_cols):
            for row in xrange(self.n_elec_rows):
                z_pos = -0.019 +col*self.electrode_pitch/np.sqrt(2)
                y_pos = 0.88+ (self.n_elec_rows -row)*self.electrode_pitch*np.sqrt(2)
                if col % 2 == 1:
                    y_pos += self.electrode_pitch/np.sqrt(2)
                self.elec_z[i], self.elec_y[i]= z_pos, y_pos
                i+=1

class HD_MEA_CMOS128_Fig:
    """HD_MEA_CMOS128 class. High-density Multielectrode
    Array with highest posible density. Honeycomb patters"""
    def __init__(self, elec_x, cell, comp):
        self.pattern = 'honeycomp'
        self.electrode_pitch = 0.017 * 8
        self.electrode_separation = self.electrode_pitch*2 #0.2/100
        self.n_elec_rows = 7
        self.n_elec_cols = 3
        self.n_elecs = self.n_elec_cols * self.n_elec_rows
        self.cylinder_height = 3.
        self.sigma_E = .3
        self.sigma_S = 3.
        self.elec_x = elec_x
        self.elec_z = np.zeros(self.n_elecs)
        self.elec_y = np.zeros(self.n_elecs)
        i = 0
        for col in xrange(self.n_elec_cols):
            for row in xrange(self.n_elec_rows):
                z_pos = -0.115 +col*self.electrode_pitch/np.sqrt(2)
                y_pos = -0.3+ (self.n_elec_rows -row)*self.electrode_pitch*np.sqrt(2)
                if col % 2 == 1:
                    y_pos += self.electrode_pitch/np.sqrt(2)
                self.elec_z[i], self.elec_y[i]= z_pos, y_pos
                i+=1



class HD_MEA_CMOS128_HAY:
    """HD_MEA_CMOS128 class. High-density Multielectrode Array with highest posible density. Honeycomb patters"""
    def __init__(self, elec_x, cell, comp):
        self.pattern = 'honeycomp'
        self.electrode_pitch = 0.017
        self.electrode_separation = self.electrode_pitch*2 #0.2/100
        self.n_elec_rows = 1
        self.n_elec_cols = 1
        self.n_elecs = self.n_elec_cols * self.n_elec_rows
        self.cylinder_height = 3.
        self.sigma_E = .3
        self.sigma_S = 3.
        self.elec_x = elec_x
        self.elec_z = np.array([cell.zmid[comp]/1000])#np.zeros(self.n_elecs)
        self.elec_y = np.array([cell.ymid[comp]/1000])#np.zeros(self.n_elecs)
        i = 0
        #for col in xrange(self.n_elec_cols):
        #    for row in xrange(self.n_elec_rows):
        #        z_pos = -0.019 +col*self.electrode_pitch/np.sqrt(2)
        #        y_pos = 0.88+ (self.n_elec_rows -row)*self.electrode_pitch*np.sqrt(2)
        #        if col % 2 == 1:
        #            y_pos += self.electrode_pitch/np.sqrt(2)
        #        self.elec_z[i], self.elec_y[i]= z_pos, y_pos
        #        i+=1


                
class SalamanderSimulationFelix:
    def __init__(self):
        self.timers_NEURON = 2**-7#2**-6
        self.timers_python = 2**-5
        self.samp_rate = 32000
        self.tstartms = 0
        self.tstopms = 10100
        self.ntsteps_org = (self.tstopms - self.tstartms)*1/self.timers_python +1
        self.ntsteps_downsampled = (self.tstopms - self.tstartms)\
                                   * self.samp_rate/1000 +1 
        self.n_box_cols = 16
        self.n_box_rows = 8
        self.n_neurons = self.n_box_cols * self.n_box_rows
        self.box_side_length = 27.
        self.box_margin = self.box_side_length/10

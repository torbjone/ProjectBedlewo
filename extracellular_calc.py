import MoI
import MEA

import numpy as np
from sys import stdout
import os
from plot_mea import plot_interval_on_all_elecs, plot_neuron_from_side,\
     animate_MEA, plot_electrodes_and_neuron
from tools import analyze_neuron

def make_mapping(cell, Mea, set_up_parameters, output_folder, do_calculation):
    if do_calculation:
        print '\033[1;35mMakeing mapping ...\033[1;m'
        t = set_up_parameters['slice_thickness']
        sigma_1 = set_up_parameters['sigma_1']
        sigma_2 = set_up_parameters['sigma_2']
        sigma_3 = set_up_parameters['sigma_3']
        comp_coors = np.array([cell.xmid/1000, cell.ymid/1000, cell.zmid/1000])
        currents =  cell.imem
        n_compartments = len(currents[:,0])
        n_elecs = Mea.n_elecs
        mapping = np.zeros((Mea.n_elecs,n_compartments))
        steps = set_up_parameters['steps']
        print "Number of compartments: %i" % n_compartments
        elec_x =  -t/2
        elec_y = Mea.elec_y
        elec_z = Mea.elec_z
        for comp in xrange(n_compartments):
            percentage = (comp+1)*100/n_compartments
            stdout.write("\r%d %% complete" % percentage)
            stdout.flush()
            #if (cell.ymid[comp] < 800) or cell.zmid[comp] > 0.00:
            #        continue
            for elec in xrange(n_elecs):
                elec_pos = [elec_x, elec_y[elec], elec_z[elec]]
                #charge_x = comp_coors[0,comp]
                charge_pos = comp_coors[:,comp]
                #print elec_pos
                #print charge_pos
                mapping[elec, comp] += Moi.anisotropic_moi(\
                    charge_pos, elec_pos)
        print ''
        np.save(output_folder + 'mapping.npy', mapping)
    else:
        mapping = np.load(output_folder + 'mapping.npy')
    return mapping

def find_signal_at_electrodes(cell, Mea, mapping, output_folder, do_calculation):
    if do_calculation:
        print '\033[1;35mFinding signal at electrodes ...\033[1;m'
        n_t_steps = len(cell.tvec)
        meaArray = np.zeros((Mea.n_elecs, n_t_steps))
        current = cell.imem
        n_compartments = len(current[:,0])
        for elec in xrange(Mea.n_elecs):
            for comp in xrange(n_compartments):
                meaArray[elec,:] += mapping[elec,comp] * current[comp, :]
        np.save(output_folder + 'signal.npy', meaArray)
    else:
        meaArray = np.load(output_folder + 'signal.npy')
    return meaArray

def make_test_cell():
    class Cell:
        def __init__(self):
            self.imem = np.array([[0,-1,1,0], [0,-1,1,0]])
            self.xmid = np.array([0, -50.])
            self.ymid = np.array([0, 700.])
            self.zmid = np.array([0, 100.])
            self.tvec = np.arange(len(self.imem[0,:]))
            self.timeres_python = 0.1
    cell = Cell()
    return cell

def make_test_mea():
    class MEA:
        def __init__(self):
            self.electrode_pitch = 0.017
            self.electrode_separation = self.electrode_pitch*2 #0.2/100
            self.n_elec_rows = 16
            self.n_elec_cols = 8
            self.n_elecs = self.n_elec_rows*self.n_elec_cols
            self.elec_z = np.zeros(self.n_elecs)
            self.elec_y = np.zeros(self.n_elecs)
            i = 0
            for col in xrange(self.n_elec_cols):
                for row in xrange(self.n_elec_rows):
                    z_pos = col*self.electrode_pitch/np.sqrt(2)
                    y_pos = (self.n_elec_rows -row)*self.electrode_pitch*np.sqrt(2)
                    if col % 2 == 1:
                        y_pos += self.electrode_pitch/np.sqrt(2)
                    self.elec_z[i], self.elec_y[i]= z_pos, y_pos
                    i+=1
    Mea = MEA()
    return Mea
set_up_parameters = {
    'sigma_1': [0.0, 0.0, 0.0], # Below Electrode
    'sigma_2': [0.3, 0.3, 0.3], # Tissue
    'sigma_3': [3.0, 3.0, 3.0], # Saline
    'slice_thickness': 0.2,
    'steps' : 20,
    'detection_limit': 30, # [mV] spike sorting limit
    }

if __name__ == '__main__':
    neuron_folder = 'larkum_tuft_spike/'
    output_folder = 'extracellular_test/'
    elec_x = -set_up_parameters['slice_thickness']/2
    Mea = MEA.HD_MEA_CMOS128(elec_x)
    #Mea = make_test_mea()
    from larkum_sim import get_cell
    cell = get_cell(neuron_folder, False)
    #cell = make_test_cell()
    Moi = MoI.MoI(set_up_parameters, True)
    #plot_electrodes_and_neuron(cell, Mea)
    mapping = make_mapping(cell, Mea, set_up_parameters, output_folder, False)
    signal = find_signal_at_electrodes(cell, Mea, mapping, output_folder, False)
    #print mapping
    #plot_interval_on_all_elecs(cell, signal, [20,30], Mea)
    #plot_neuron_from_side(cell, Mea, set_up_parameters)
    animate_MEA(signal, cell, Mea, [5,10])
    #analyze_neuron(cell, mapping, Mea, signal_range = [5,10])


import matplotlib.pylab as pl
import numpy as np
import LFPy
import os
import sys
import neuron
#seed for random generation
pl.seed(0)

def push_simulation_to_folder(save_to_folder, data_from_folder):
    print "Copying all simulation results from %s " %data_from_folder\
          +"and simulation file to folder %s." % save_to_folder
    try:
        os.mkdir(save_to_folder)
    except(OSError):
        print "Result folder already exists. Overwriting..."
    os.system('cp %s/* %s' %(data_from_folder, save_to_folder))


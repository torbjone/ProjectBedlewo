#!/usr/bin/env python
"""minimal mods to Espen's script to generate polytrode data"""
import pylab as pl
import numpy as np
import nest
nest.set_verbosity('M_FATAL')
from mpi4py import MPI
import os, sys
import LFPy
import ipdb

# add testdata
try:
    import testdata
except:
    home = os.path.expanduser('~')
    path = os.path.join(home, 'umb', 'CRCNS/trunk')
    if not path in sys.path: sys.path.append(path)

from testdata.cellsimmethods import \
    cellsim_active, draw_rand_pos, shufflemorphos, \
    shufflecustom_codes, collect_data


from testdata.nestfun import run_brunel_delta_nest, gdfFilesProcessing
from time import time

#set some seeds
SEED = 12345678
NESTSEED = SEED
pl.seed(SEED)

################# Initialization of MPI stuff ##################################
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()
MASTER_MODE = COMM.rank == 0
print 'SIZE %i, RANK %i, MASTER_MODE: %s' % (SIZE, RANK, str(MASTER_MODE))


#print out memory consumption etc every ten seconds
if RANK == 0 or RANK == 8:
    if sys.platform == 'darwin':
        pass
    else:
        os.system("vmstat 10 -S M &")


################################################################################
# PARAMETERS
################################################################################
cellParams = {
    'v_init' : -80,
    'passive' : False,
    'nsegs_method' : 'lambda_f',
    'lambda_f' : 100,
    'timeres_NEURON' : 2**-5,
    'timeres_python' : 2**-5,
    'tstartms' : -500, #-500.,
    'tstopms' : 1000, # 10000.
    'verbose' : True,
}
simulationParams = {
}

POPULATION_SIZE = 20
# populationParams = {
#     'number' : POPULATION_SIZE,
#     'radius' : 50,
#     'z_min' : -75,
#     'z_max' : 25
# }

cellpositions = []
pos = np.zeros((POPULATION_SIZE, 3))
for i in range(POPULATION_SIZE):
    z = (np.random.rand()-.5)*2000
    x, y = np.random.randn(2)*10
    cellpositions.append({'xpos' : x, 'ypos' : y, 'zpos' : z})
    pos[i] = x,y,z

ipdb.set_trace()

# electrodeParams = {
#     'x' : pl.array([0,-1,pl.sin(pl.pi/6),pl.sin(pl.pi/6)])*25.,
#     'y' : pl.array([0,0,-pl.cos(pl.pi/6),pl.cos(pl.pi/6)])*25.,
#     'z' : pl.array([-50.,0,0,0]),
#     'sigma' : 0.1,
#     'N' : pl.array([   [0,0,-1],[-1*pl.cos(pl.pi/9),0,-1*pl.sin(pl.pi/9)],
#                     [pl.sin(pl.pi/6)*pl.cos(pl.pi/9),
#                         -pl.cos(pl.pi/9)*pl.cos(pl.pi/9),-1*pl.sin(pl.pi/9)],
#                     [-pl.sin(pl.pi/6)*pl.cos(pl.pi/9),
#                         -pl.cos(pl.pi/9)*pl.cos(pl.pi/9),1*pl.sin(pl.pi/9)]]),
#     'r' : 7.,
#     'n' : 100,
#     'r_z': pl.array([[-1E199,-50.00001,-50,75,1E99],[0,0,7,48,48]]),
#     'seedvalue' : None,
# }

ch = 32
N = pl.empty((ch, 3))
for i in xrange(N.shape[0]): N[i,] = [1, 0, 0] #normal unit vec. to contacts
electrodeParams = {             #parameters for electrode class
    'sigma' : 0.1,              #Extracellular potential
    'x' : pl.zeros(ch),         #Coordinates of electrode contacts
    'y' : pl.zeros(ch),
    'z' : pl.linspace(-500,1000,ch),
    'n' : 20,
    'r' : 15,
    #'r_z' : pl.array([[-1E199,-500.00001,-500.00001,75,1E99],[0,0,7,48,48]]),
    'N' : N,
}

synparams_AMPA = {         #Excitatory synapse parameters
    'e' : 0,           #reversal potential
    'syntype' : 'Exp2Syn',   #conductance based exponential synapse
    'tau1' : 1.,         #Time constant, rise
    'tau2' : 3.,         #Time constant, decay
    'weight' : 0.0020,   #Synaptic weight
    'color' : 'r',      #for pl.plot
    'marker' : '.',     #for pl.plot
}
synparams_GABA_A = {         #Inhibitory synapse parameters
    'e' : -80,
    'syntype' : 'Exp2Syn',
    'tau1' : 1.,
    'tau2' : 12.,
    'weight' : 0.0020,
    'color' : 'b',
    'marker' : '.',
}

#where to insert, densitiy of synapses
insert_synapses_AMPA_args = {
    'section' : 'alldend',
    'nPerArea' : 12.5E-3,
}
insert_synapses_GABA_A_args = {
    'section' : 'alldend',
    'nPerArea' : 10.0E-3,
}

#parameter for NEST simulation
nestSimParams = {
    'simtime' :     cellParams['tstopms']-cellParams['tstartms'],
    'dt' :          cellParams['timeres_NEURON'],
    'order' :       2500,
    'g' :           5.0,
    'eta' :         2.0,
    'epsilon' :     0.05,
    'J' :           0.1,
    #'grng_seed' :   NESTSEED,
}


noiseParams = {
    'T' : cellParams['tstopms'] * 1E-3,
    'dt' : cellParams['timeres_python'] * 1E-3,
    'sigma' : 1.0,
} 
NOISECORRELATION = 0.9

#nyquist frequency
nyquist = 1000. / cellParams['timeres_python'] / 2 

#template length and offset
TEMPLATELEN = 200
OFFS = 0.3
    
# set the rotations
rotations = []
defaultrotation = {}

savelist = [
    'AP_train',
    'somav',
    'timeres_NEURON',
    'timeres_python',
    'somapos',
    'x',
    'y',
    'z',
    'LFP',
    'morphology',
    'default_rotation',
    'custom_code',
    'electrodecoeff',
]

morphologies = [
    'morphologies/cell1.hoc',
    'morphologies/cell2.hoc',
    'morphologies/cell3.hoc',
]

custom_codes = [
    ['custom_codes.hoc', 'biophys1.hoc'],
    ['custom_codes.hoc', 'biophys2.hoc'],
    ['custom_codes.hoc', 'biophys3.hoc'],
    #['custom_codes.hoc', 'biophys4.hoc'],
]

################################################################################
## MAIN
################################################################################

TIME = time()

##### load mechanisms ##########################################################
LFPy.cell.neuron.load_mechanisms("../mod")

##### NEST SIMULATION FOR SPIKE TIMES, RUN IN PARALLELL ACROSS RANKS ##########
#delete old files from prev nest sim
if MASTER_MODE:
    os.system('rm -r *.gdf')

run_brunel_delta_nest(**nestSimParams)
ipdb.set_trace()

#process gdf files, creating savedata/SpTimes*.h5 files with spike times
gdfFilesProcessing(nestSimParams)


#extract so many spikes for each cell for excitatory and inhib nest-cells
if MASTER_MODE:
    NTIMES = 1000
    SpCell0 = pl.randint(0, nestSimParams['order']*4,
                         size=(POPULATION_SIZE, NTIMES)).astype('int32')
    SpCell1 = pl.randint(0, nestSimParams['order'],
                         size=(POPULATION_SIZE, NTIMES)).astype('int32')
else:
    SpCell0 = None
    SpCell1 = None

#broadcast the cell indexes to all ranks
SpCell0 = COMM.bcast(SpCell0, root=0)
SpCell1 = COMM.bcast(SpCell1, root=0)

print 'created SpTimes in %.6f' % (time() - TIME)

#distrinite soma locations and the shuffled morphos etc
if MASTER_MODE:
    ##### create lists/dicts of locations, morphologies, codes #################
    # pop_soma_pos = draw_rand_pos(min_r = electrodeParams['r_z'],
    #                              **populationParams)
    pop_soma_pos = cellpositions
    shuffledmorphologies = shufflemorphos(morphologies)
    shuffledcustom_codes = shufflecustom_codes(custom_codes)
    for cellindex in xrange(POPULATION_SIZE):
        defaultrotation.update({'z' : pl.random() * 2 * pl.pi})
        rotations.append(defaultrotation.copy())
else:
    pop_soma_pos = None
    shuffledmorphologies = None
    shuffledcustom_codes = None
    rotations = None

#broadcasting
pop_soma_pos = COMM.bcast(pop_soma_pos, root=0)
shuffledmorphologies = COMM.bcast(shuffledmorphologies, root=0)
shuffledcustom_codes = COMM.bcast(shuffledcustom_codes, root=0)
rotations = COMM.bcast(rotations, root=0)

################################################################################
# WORKERS, MASTER
################################################################################

#initialize workers
if not MASTER_MODE and SIZE > 1:
    while(True):
        #receive parameters from master
        cellindex = COMM.recv(source=0)
        if cellindex == None:
            break
        else:
            #send simulation results to master
            COMM.send(cellsim_active(cellindex, cellParams, simulationParams,
                            electrodeParams,
                            shuffledmorphologies, shuffledcustom_codes,
                            synparams_AMPA, synparams_GABA_A,
                            insert_synapses_AMPA_args,
                            insert_synapses_GABA_A_args,
                            SpCell0[cellindex], SpCell1[cellindex],
                            nestSimParams,
                            pop_soma_pos, rotations, savelist), dest=0)


#master will send parameters to workers and receive simulation results unordered
if MASTER_MODE:
    #init container for simulation results
    outputs = []
    #using MPI
    if SIZE > 1:
        dest = 1 #counter on [1, SIZE]
        #sending parameters to workers
        for cellindex in xrange(POPULATION_SIZE):
            COMM.send(cellindex, dest=dest)
            dest += 1
            #reset counter:
            if dest >= SIZE:
                dest = 1

        #receiving simulation results from any worker, storing in container
        for cellindex in xrange(POPULATION_SIZE):
            COMM.recv(source=MPI.ANY_SOURCE)
            
        #killing workers
        for i in xrange(1, SIZE):
            COMM.send(None, dest=i)
    
    #serial mode
    else:
        for cellindex in xrange(POPULATION_SIZE):
            cellsim_active(cellindex, cellParams,
                simulationParams, electrodeParams,
                shuffledmorphologies, shuffledcustom_codes,
                synparams_AMPA, synparams_GABA_A,
                insert_synapses_AMPA_args, insert_synapses_GABA_A_args,
                SpCell0[cellindex], SpCell1[cellindex],
                nestSimParams,
                pop_soma_pos, rotations, savelist)
            
    print 'execution time: %.3f seconds' %  (time()-TIME)


################################################################################
# COLLECT DATA FOR EACH CELL ON MASTER_
################################################################################

if MASTER_MODE:
    populationParams = None
    cells = collect_data(cellParams, POPULATION_SIZE,
                 shuffledmorphologies,
                 populationParams, pop_soma_pos, rotations,
                 noiseParams, NOISECORRELATION,
                 electrodeParams,
                 TEMPLATELEN, OFFS, nyquist)
    
    print 'Simulation done in %s seconds' % (time() - TIME)

#print out memory consumption etc every ten seconds
if RANK == 0 or RANK == 8:
    if sys.platform == 'darwin':
        pass
    else:
        os.system("killall vmstat")

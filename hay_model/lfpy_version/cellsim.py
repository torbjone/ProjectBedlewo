#!/usr/bin/env python
import pylab as pl
import LFPy

#load compiled mechs from the mod-folder
LFPy.cell.neuron.load_mechanisms("../mod")


pl.interactive(1)



cellParams = {
    'morphology' : 'morphologies/cell1.hoc',
    #'morphology' : '../morphologies/cell1.asc',
    'passive' : False,
    'custom_code' : ['custom_codes.hoc', 'biophys3_modified.hoc'],
    'nsegs_method' : None,
    'v_init' : -80,
    'tstartms' : 0,
    'tstopms' : 3000,
    'timeres_NEURON' : 2**-3,
    'timeres_python' : 2**-3,
    'verbose' : False,
}

cell = LFPy.Cell(**cellParams)

PointProcParams = {
    'idx' : 0,
    'record_current' : False,
    'pptype' : 'IClamp',
    'amp' : 0.793,
    'dur' : 2000,
    'delay' : 700,
}

pointProcess = LFPy.StimIntElectrode(cell, **PointProcParams)

cell.simulate(rec_variables = [])

pl.plot(cell.tvec, cell.somav)
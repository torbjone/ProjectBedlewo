Author: Etay Hay, 2011
  Models of Neocortical Layer 5b Pyramidal Cells Capturing a Wide Range of
  Dendritic and Perisomatic Active Properties
  (Hay et al., PLoS Computational Biology, 2011) 

NEURON models and model sets corresponding to the paper:

folder: models
==============
NEURON code models, shown in various figures

L5PCbiophys1 - figure 1 (constrained only for BAC firing)
L5PCbiophys2 - figure 2 (constrained only for current step firing)
L5PCbiophys3 - figure 4 (constrained both for BAC and current step firing)
L5PCbiophys4 - figure S5 (AP initiation at the axon)
L5PCtemplate - general cell template

folder: simulation code
=======================
simulation code for BAC firing or step current firing.

folder: model sets
==================
model sets corresponding to various figures.
 
models_errors file: error values matrix (rows: models; columns: objectives)
models_parameters file: parameter values matrix (rows: models; columns: parameters)
objetives file: objective names, each corresponding to a column in models_errors file
genome file: parameter names and search limits, each corresponding to a column in
   the models_parameters file
   
folder: mechanisms
==================
mod files of the conductance mechanisms used.

folder: morphologies
====================
The three morphologies used in the paper. 



# SCAT
Spectrogram correlation and transformation model inspired by echolocating bats' auditory system
# Directory Contents
files
* sandbox.m -- add the subfolders into the current path of MATLAB
* freqHopping.m -- the main program to run 

subfolders:
 - coch_filters - filterbank functions to model basilar membrane motion
 - config_files - Supplementary configuration data files
 - neur_discharge - Acoustic-to-neural transduction models by IHCs
 - neur_temporal - Temporal processing models and delay lines 
 - exampleSignals - Signals that are processed by the model
 - matlib_lib - Additional libraries
 - newfuncs - functions created after initial version from Jason
 - wavParameters - Contains the wavParam structures that will be used in function: linear_...
 
 This is an extension of code from Jason Gaudette (link: https://github.com/gaudetteje/biscat). Please also see more references and short development history of this model from the above link
 
 # How to run the code?
 First run sandbox.m to add all necessary subfolders to the current MATLAB path, then run freqHopping.m
 # Contact
 Chen Ming (chen_ming@brown.edu)

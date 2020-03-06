# SCAT
Spectrogram correlation and transformation (SCAT) model, an efficient signal processing method inspired by echolocating bats' auditory system
# Directory Contents
## files
* sandbox.m -- add the subfolders into the current path of MATLAB
* freqHopping.m -- run this program to demonstrate frequency hopping.
* scatModelDemo.m -- the main program to run through the SCAT model, outputs three figures. 1. the output of filterbank; 2. the dechirped image; 3. the triangular network

## subfolders:
 - coch_filters - filterbank functions to model basilar membrane motion
 - config_files - Supplementary configuration data files
 - neur_discharge - Acoustic-to-neural transduction models by IHCs
 - neur_temporal - Temporal processing models and delay lines 
 - exampleSignals - Signals that are processed by the model
 - matlib_lib - Additional libraries
 - newfuncs - functions created after initial version from Jason
 - wavParameters - Contains the wavParam structures that will be used in function: linear_...
 
## Highlight of SCAT model: amplitude latency trading (ALT) effect
 SCAT model mimicks the auditory system of big brown bats. By adding ALT effect, the notches in dechirped image are magnified and easier to be found.
 
![](/figures/DechirpedImage-1.png)

![](/figures/DechirpedImage-2.png)
 
 This is an extension of code from Jason Gaudette (link: https://github.com/gaudetteje/biscat). Please also see more references and short development history of this model from the above link
 
 # How to run the code?
 First run sandbox.m to add all necessary subfolders to the current MATLAB path
 * then run freqHopping.m to demonstrate the part of SCAT model that can exclude the "phantom" echoes
 * run scatModelDemo.m to demonstrate the main part of SCAT model
 # Contact
 Chen Ming (chen_ming@brown.edu)


# SCAT
Spectrogram correlation and transformation (SCAT) model, an efficient signal processing method inspired by echolocating bats' auditory system. SCAT model runs a pulse-echo pair through a gammatone filterbank to get the time-frequency representations of the signal in parallel frequency channels as the first step.

![](/figures/filterbankOutput.png)

Then SCAT analyses each channel to get the time delay of target, and the geometry information of the target, i.e., fine delay between glints (reflecting points) in the target. By adding up information from all frequency channels, SCAT model can eventually make the best estimate of the two delays mentioned above.

# Directory Contents

## PLoS Comp. Biol. submission:
* bat2HFM.m -- the main code to demonstrate the module in SCAT model for processing bat 2-harmonic FM sweeps.
* clutterRejection.m -- the code to demo the ALT could further delay the arrival of off-axis clutter echoes.
* dolClick.m -- SCAT model can also be used to process dolphin clicks.
* clutterRejection_dolClick.m -- Clutter rejection example for dolphin clicks

Each of the four programs will take a few minutes to run.

## files
* sandbox.m -- add the subfolders into the current path of MATLAB
* freqHopping.m -- run this program to demonstrate frequency hopping.
* scatModelDemo.m -- the main program to run through the SCAT model, outputs three figures. 1. the output of filterbank; 2. the dechirped image; 3. the triangular network

*(updated on 6/12/2020)*:
* bat2HFM.m -- the main code to demonstrate the module in SCAT model for processing bat 2-harmonic FM sweeps.
* clutterRejection.m -- the code to demo the ALT could further delay the arrival of off-axis clutter echoes.
* dolClick.m -- SCAT model can also be used to process dolphin clicks.
* binauralTracking.m -- An addition of binaural tracking was added based on SCAT model. It's a demo for how a active sonar model can find the right target.

*(updated on 6/23/2020)*:
* clutterRejection_dolClick.m -- Clutter rejection example for dolphin clicks


## subfolders:
 - coch_filters - filterbank functions to model basilar membrane motion
 - config_files - Supplementary configuration data files
 - neur_discharge - Acoustic-to-neural transduction models by IHCs
 - neur_temporal - Temporal processing models and delay lines 
 - exampleSignals - Signals that are processed by the model
 - matlib_lib - Additional libraries
 - newfuncs - functions created after initial version from Jason
 - wavParameters - Contains the wavParam structures that will be used in function: linear_...
 

 
# Highlight of SCAT model: amplitude latency trading (ALT) effect
 SCAT model mimicks the auditory system of big brown bats. By adding ALT effect, the notches in dechirped echo are magnified and easier to be found.
 
![](/figures/DechirpedImage-1.png)

![](/figures/DechirpedImage-2.png)
 
 This is an extension of code from Jason Gaudette (link: https://github.com/gaudetteje/biscat). Please also see more references and short development history of this model from the above link
 
 # How to run the code?
 First run sandbox.m to add all necessary subfolders to the current MATLAB path
 * then run freqHopping.m to demonstrate the part of SCAT model that can exclude the "phantom" echoes
 * run scatModelDemo.m to demonstrate the main part of SCAT model
 # Contact
 Chen Ming (chen_ming@brown.edu)
 
 I am a postdoc at Neuroscience department, Brown University. I work with Prof. James Simmons. Please feel free to reach out if you have questions or suggestions.
 
 I am writing a journal paper that describes this model. Please check back later to get a link of the paper. Also, I am constantly updating the model. Please check back and download the latest one.


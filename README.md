# DrivingFatigueDetection

This is a Matlab implementation of our paper "Sequential Nonlinear Encoding; A Low Dimensional Regression Algorithm with Application to EEG based Driving Fatigue Detection"
# Prerequisite
Matlab R2017a, Windows 7 (64 bits).\
Wavelet toolbox.\
Signal toolbox.\
Curefit toolbox.
# Data preparation
     We use the dataset provided by WL Zheng (Center for Brain-like Computing and Machine Intelligence, Shanghai Jiao Tong
     University), including EEG samples of 23 subjects recorded by a Neuroscan system at a sampling rate of 1000 Hz. Each
     sample includes 17 channels acquired from the posterior (CP1, CPZ, CP2, P1, PZ, P2, PO3, POZ, PO4, O1, OZ, and O2) and
     temporal (FT7, FT8, T7, T8, TP7, and TP8) areas of brain according to the international 10–20 electrode system.

All the processing tasks including artifact removal, osilatory pattern separation, epoch dissection, and logarithmic feature extraction are aggregated into a single operation "featureextraction.m".

For evaluation, please run the featureextraction.m and save the output with a .mat extension.

Then, apply the regression task using main.m which simultanously outputs the regression curves for our training and test sessions.
 
![alt text](https://github.com/94147852369/DrivingFatigueDetection/blob/master/topo.png)

# Acknowledgement

Our implementation borrowed some functions from original Matlab toolboxes.


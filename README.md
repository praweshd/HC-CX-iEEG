# HC-CX-iEEG
EEG Analysis Toolbox
This repository contains a set of MATLAB scripts for the analysis of EEG data. The scripts cover various aspects of EEG signal processing, event detection, sleep scoring, and data cleaning. Below is an overview of the functionalities provided by each script.

1. EEG Stream and Downsampling (dat2lfp.m)
Streams EEG data sampled at 20 kHz in 60-second durations and performs downsampling to 1250 Hz.

2. Data Tracker & Loader (Dat_tracker.m)
Loads a specific segment of EEG data centered around a specified time for a certain duration.

3. IED (25-80 Hz) Pathological Discharges Detection (IED_Detect_mof_v2.m)
Detects Interictal Epileptiform Discharges (IEDs) on all electrodes using a combination of frequency, amplitude, and duration parameters. The process includes bandpass filtering, signal rectification, and event elimination based on specific criteria.

4. Step-by-step Signal Processing for EEG Event Detection - Tester Code (gamma_tester.m)
Provides a step-by-step demonstration of signal processing for EEG event detection.

5. Spindle Detections (Spindle_Detect_wavelet.m + FindSpindles_wav.m)
Detects spindle events based on wavelet-derived power and duration parameters. The approach involves analyzing the ratio of normalized autoregressive wavelet-based power in different frequency bands.

6. Sleep Scoring (SleepScoreMaster.m + ClusterState_GetMetrics.m)
Resamples electrophysiological data to 1250 Hz for local field potential (LFP) analysis. Identifies sleep epochs based on accelerometer motion signals, absence of electromyogram artifacts, and specific frequency characteristics in the intracranial electroencephalography (iEEG) spectrograms.

7. Cleaning NREM (CleanNREM.m)
Cleans artifacts within the NREM sleep section post-sleep scoring by assigning them as WAKE.

8. LFP Matching (LFPMatching.m)
Identifies and isolates IED events in channel 2 (res2) that are matched/unmatched with events in channel 1 (res1) within a certain tolerance (time in ms).

9. Hippocampal-Ripple Triggered Power in Cortical Channels (HippRipplePower.m)
Analyzes the power of cortical channels triggered by hippocampal ripple events.

10. Identifying Animal-Head Position from LED Tracking Video (LEDPositionIdentification.m)
Uses LED tracking video to identify the position of the animal's head.

11. Perform PCA/K-means for SPI Group Clustering (PCA_KMeans_SPI.m)
Applies Principal Component Analysis (PCA) and K-means clustering to investigate the presence of two clusters of SPI (Spindle) groups based on local and global power.

Feel free to explore and utilize these scripts for EEG data analysis. If you have any questions or suggestions, please don't hesitate to reach out!

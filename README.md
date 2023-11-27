# HC-CX-iEEG

EEG Analysis Toolbox
This repository contains a set of MATLAB scripts for the analysis of EEG data. The scripts cover various aspects of EEG signal processing, event detection, sleep scoring, and data cleaning. Below is an overview of the functionalities provided by each script.

## 1. **EEG Stream and Downsampling (`dat2lfp.m`)**
Streams EEG data sampled at 20 kHz in 60-second durations and performs downsampling to 1250 Hz.

## 2. **Data Tracker & Loader (Dat_tracker.m)**
Loads a specific segment of EEG data centered around a specified time for a certain duration.

## 3. **IED (25-80 Hz) Pathological Discharges Detection (IED_Detect_mof_v2.m)**
Detects Interictal Epileptiform Discharges (IEDs) on all electrodes using a combination of frequency, amplitude, and duration parameters: (i) bandpass filtering at 25–80 Hz and signal rectification; (ii) detection of events where the filtered envelope was >3 times above baseline; (iii) elimination of events where the unfiltered envelope was <3 times above baseline; and (iv) elimination of IEDs occurring within 500 ms of another IED to prevent over-correlation due to a run of IEDs


## 4. **Step-by-step Signal Processing for EEG Event Detection - Tester Code (gamma_tester.m)**
Provides a step-by-step demonstration of signal processing for EEG event detection.

## 5. **Spindle Detections (Spindle_Detect_wavelet.m + FindSpindles_wav.m)**
Detects spindle events based on wavelet-derived power and duration parameters. A ratio of normalized autoregressive wavelet-based P_AR where spindle band power (P_spi) was based on 10–20 Hz, low band power was based on 2–8 Hz (P_low), and high band power was based on 25–40 Hz (P_high). Spindle events were identified when the ratio crossed zero and was >0.1 for a minimum of 300 ms and a maximum of 3 s. We detected both fast (13–15 Hz) and slow (9–12 Hz) spindles using this approach. All detections were visually inspected for accuracy for each recording session.

## 6. **Sleep Scoring (SleepScoreMaster.m + ClusterState_GetMetrics.m)**
Identifies sleep epochs based on accelerometer motion signals, absence of electromyogram artifacts, and specific frequency characteristics in the intracranial electroencephalography (iEEG) spectrograms.

The electrophysiological data were resampled to 1250 Hz to facilitate local field potential (LFP) analysis. We identified epochs of sleep by first detecting immobility in the motion signal of the accelerometer attached to the headstage and absence of electromyogram (EMG) artifacts. NREM sleep epochs were then detected by locating periods of elevated delta (0.5 to 4 Hz) amplitude in the neocortex at times of immobility. rapid eye movement sleep (REM) epochs were distinguished by an increased theta–delta-band frequency ratio in the intracranial electroencephalography (iEEG) spectrograms in the hippocampus. Sleep-scored epochs were then visually inspected and manually adjusted using whitened spectrograms and raw traces to eliminate short epochs containing movement artifacts.


## 7. **Cleaning NREM (CleanNREM.m)**
Cleans artifacts within the NREM sleep section post-sleep scoring by assigning them as WAKE.

## 8. **LFP Matching (LFPMatching.m)**
Identifies and isolates IED events in channel 2 (res2) that are matched/unmatched with events in channel 1 (res1) within a certain tolerance (time in ms).

## 9. **Hippocampal-Ripple Triggered Power in Cortical Channels**
Analyzes the power of cortical channels triggered by hippocampal ripple events.

## 10. **Identifying Animal-Head Position from LED Tracking Video** 
Uses LED tracking video to identify the position of the animal's head.

## 11. **Perform PCA/K-means for SPI Group Clustering**  
Applies Principal Component Analysis (PCA) and K-means clustering to investigate the presence of two clusters of SPI (Spindle) groups based on local and global power.

 

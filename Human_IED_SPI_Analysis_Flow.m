% Human IED-SPI Analysis Flow

%% Load File Information

CH_key=dir('*CH_key.mat');
load (CH_key.name)
CH_N=length (CH_key);

lfp_file=dir('*.lfp');
lfp_filename=lfp_file.name;

%% IED Detection

for ii = 1:length(CH_key)
tic
IED_Detect_mod(CH_N, ii,'med_thresh',3, 'mean_thresh',3);
toc
end

%% Spindle Detection

tic
for SPI_CH = 1:length(CH_key)         

    Spindle_Detect(CH_N,SPI_CH);    % thresholds [2 4]
    
end
toc

%% Eliminate IEDs from Spindle Detection

RemoveIED_Spidetection          % Removes all spindles that have an IED detected on the same CH within the start and end times

%% Create IED, SPI, and SPIONLY Structures

IED_SPImat2struct

%% 


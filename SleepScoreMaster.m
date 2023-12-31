function SleepState = SleepScoreMaster(basePath,EMGCH,ACCM,varargin)
%SleepScoreMaster(datasetfolder,recordingname)
%This is the master function for sleep state scoring.
%
%INPUT 
%   basePath        folder containing .xml and .lfp files.
%                   basePath and files should be of the form:
%                   'whateverfolder/recordingName/recordingName'
%   (optional)      If no inputs included, select folder(s) containing .lfp
%                   and .xml file in prompt.
%   (optional)      if no .lfp in basePath, option to select multiple 
%                   lfp-containing subfolders
%                          
%   OPTIONS
%   'savedir'       Default: datasetfolder
%   'overwrite'     Default: false, overwrite all processing steps
%   'savebool'      Default: true
%   'scoretime'     Default: [0 Inf] NOTE: must be continous interval until
%                   someone updates this...
%   'SWWeightsName' Name of file in path (in Dependencies folder) 
%                   containing the weights for the various frequencies to
%                   be used for SWS detection.  Default is to use Power Spectrum Slope ('PSS'),
%                   but can also try 'SWweights.mat'
%                     - For hippocampus-only recordings, enter
%                     'SWweightsHPC.mat' for this
%   'Notch60Hz'     Boolean 0 or 1.  Value of 1 will notch out the 57.5-62.5 Hz
%                   band, default is 0, no notch.  This can be necessary if
%                   electrical noise.
%   'NotchUnder3Hz' Boolean 0 or 1.  Value of 1 will notch out the 0-3 Hz
%                   band, default is 0, no notch.  This can be necessary
%                   due to poor grounding and low freq movement transients
%   'NotchHVS'      Boolean 0 or 1.  Value of 1 will notch the 4-10 and 
%                   12-18 Hz bands for SW detection, default is 0, no 
%                   notch.  This can be useful in
%                   recordings with prominent high voltage spindles which
%                   have prominent ~16hz harmonics
%   'NotchTheta'    Boolean 0 or 1.  Value of 1 will notch the 4-10 Hz
%                   band for SW detection, default is 0, no notch.  This 
%                   can be useful to
%                   transform the cortical spectrum to approximately
%                   hippocampal, may also be necessary with High Voltage
%                   Spindles
%   'stickytrigger' Implements a "sticky" trigger for SW/EMG threshold 
%                   crossings: metrics must reach halfway between threshold
%                   and opposite peak to count as crossing (reduces
%                   flickering, good for HPC recordings) (default:false)
%   'SWChannels'    A vector list of channels that may be chosen for SW
%                   signal
%   'ThetaChannels' A vector list of channels that may be chosen for Theta
%                   signal
%   'rejectChannels' A vector of channels to exclude from the analysis
%   'noPrompts'     (default:false) an option to not prompt user of things
%
%OUTPUT 
%   !THIS IS OUT OF DATE - UPDATE!
%   StateIntervals  structure containing start/end times (seconds) of
%                   NREM, REM, WAKE states and episodes. states is the 
%                   "raw" state scoring. episodes are joined episodes of 
%                   extended (40s) time in a given states, allowing for 
%                   brief interruptions. also contains NREM packets, 
%                   unitary epochs of NREM as described in Watson et al 2016.
%                   saved in a .mat file:
%                   recordingname_SleepScore.mat 
%   
%
 %% Recording Selection
%if recname is 'select' or something
%use uigetfile to pick and get list of filenames
%if recname is 'all', get all recordings in a folder and
%then run SleepScoreMaster on each of the filenames'
%if no input arguements... select uigetfile

%Select from no input
if ~exist('basePath','var')
    basePath = uigetdir(cd,...
        'Which recording(s) would you like to state score?');
    if isequal(basePath,0);return;end  
end

%Separate datasetfolder and recordingname
[datasetfolder,recordingname,extension] = fileparts(basePath);
recordingname = [recordingname,extension]; % fileparts parses '.' into extension


%% If there is no .lfp in basePath, choose (multiple?) folders within basePath.
%Select from dataset folder - need to check if .xml/lfp exist
if ~exist(fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file') && ...
    ~exist(fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    disp(['no .lfp file in basePath, pick a selection of session folders',...
             'containing .lfp files'])
        %foldercontents = dir(basePath);
        %possiblerecordingnames = {foldercontents([foldercontents.isdir]==1).name};
        [basePaths,recordingname] = bz_FindBasePaths(basePath,'select',true); %Find all basePaths within the topPath
%         [s,v] = listdlg('PromptString','Which recording(s) would you like to state score?',...
%                         'ListString',baseNames);
%         recordingname = baseNames(s);
%         basePaths = basePaths(s);
        
end

%If multiple recordings, loop calling SleepScoreMaster with each
numrecs = length(recordingname);
if numrecs > 1 && iscell(recordingname)
    disp(['Multiple Recordings (',num2str(numrecs),')'])
    for rr = 1:numrecs
        multibasepath = basePaths{rr};
        SleepScoreMaster(multibasepath,varargin{:})
        close all
    end
    return
elseif numrecs == 1 && iscell(recordingname)
        recordingname = recordingname{1};
end

display(['Scoring Recording: ',recordingname]);

%% inputParse for Optional Inputs and Defaults
p = inputParser;

defaultOverwrite = false;    %Pick new and Overwrite existing ThLFP, SWLFP?
defaultSavebool = true;    %Save Stuff (EMG, LFP)

defaultSavedir = datasetfolder;

defaultScoretime = [0 Inf];
%defaultSWWeightsName = 'SWweights.mat';
defaultSWWeightsName = 'PSS';
defaultNotch60Hz = 0;
defaultNotchUnder3Hz = 0;
defaultNotchHVS = 0;
defaultNotchTheta = 0;
defaultSWChannels = 0;
defaultThetaChannels = 0;

addParameter(p,'overwrite',defaultOverwrite)
addParameter(p,'savebool',defaultSavebool,@islogical)
addParameter(p,'savedir',defaultSavedir)
addParameter(p,'scoretime',defaultScoretime)
addParameter(p,'SWWeightsName',defaultSWWeightsName)
addParameter(p,'Notch60Hz',defaultNotch60Hz)
addParameter(p,'NotchUnder3Hz',defaultNotchUnder3Hz)
addParameter(p,'NotchHVS',defaultNotchHVS)
addParameter(p,'NotchTheta',defaultNotchTheta)
addParameter(p,'SWChannels',defaultSWChannels)
addParameter(p,'ThetaChannels',defaultThetaChannels)
addParameter(p,'rejectChannels',[]);
addParameter(p,'noPrompts',false);
addParameter(p,'stickytrigger',false);

parse(p,varargin{:})
%Clean up this junk...
overwrite = p.Results.overwrite; 
savedir = p.Results.savedir;
scoretime = p.Results.scoretime;
SWWeightsName = p.Results.SWWeightsName;
Notch60Hz = p.Results.Notch60Hz;
NotchUnder3Hz = p.Results.NotchUnder3Hz;
NotchHVS = p.Results.NotchHVS;
NotchTheta = p.Results.NotchTheta;
SWChannels = p.Results.SWChannels;
ThetaChannels = p.Results.ThetaChannels;
rejectChannels = p.Results.rejectChannels;
noPrompts = p.Results.noPrompts;
stickytrigger = p.Results.stickytrigger;

%% Parameter setting
% Min Win Parameters (s): basic detection paramaters (seconds)
MinTimeWindowParms.minSWSsecs = 6;
MinTimeWindowParms.minWnexttoREMsecs = 6;
MinTimeWindowParms.minWinREMsecs = 6;       
MinTimeWindowParms.minREMinWsecs = 6;
MinTimeWindowParms.minREMsecs = 6;
MinTimeWindowParms.minWAKEsecs = 6;
%% Database File Management 
savefolder = fullfile(savedir,recordingname);
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

%Filenames of metadata and SleepState.states.mat file to save
sessionmetadatapath = fullfile(savefolder,[recordingname,'.SessionMetadata.mat']);
%Buzcode outputs
bz_sleepstatepath = fullfile(savefolder,[recordingname,'.SleepState.states.mat']);



%% Get channels not to use
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',noPrompts);
% check that SW/Theta channels exist in rec..
if length(SWChannels) > 1 
    if sum(ismember(SWChannels,sessionInfo.channels)) ~= length(SWChannels)
        error('some of the SW input channels dont exist in this recording...?')
    end   
end
if length(ThetaChannels) > 1 
    if sum(ismember(ThetaChannels,sessionInfo.channels)) ~= length(ThetaChannels)
        error('some of the theta input channels dont exist in this recording...?')
    end   
end

if exist(sessionmetadatapath,'file')%bad channels is an ascii/text file where all lines below the last blank line are assumed to each have a single entry of a number of a bad channel (base 0)
    load(sessionmetadatapath)
    rejectChannels = [rejectChannels SessionMetadata.ExtracellEphys.BadChannels];
elseif isfield(sessionInfo,'badchannels')
    rejectChannels = [rejectChannels sessionInfo.badchannels]; %get badchannels from the .xml
else
%     disp('No baseName.SessionMetadata.mat, no badchannels in your xml - so no rejected channels')
%Bypass Prawesh
    aaa=1;
end


%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
% Load/Calculate EMG based on cross-shank correlations 
% (high frequency correlation signal = high EMG).  
% Schomburg E.W. Neuron 84, 470?485. 2014)

%Prawesh - bypassing EMG giving specific channels 
%EMGCH = [160 162 175 177 173];
% EMGCH = [131,133,150,154,156]; %Jose Baseline

EMGFromLFP = bz_EMGFromLFP(basePath,EMGCH,'overwrite',overwrite,...
                                     'rejectChannels',rejectChannels,'noPrompts',noPrompts);

%% DETERMINE BEST SLOW WAVE AND THETA CHANNELS
%Determine the best channels for Slow Wave and Theta separation.
%Described in Watson et al 2016, with modifications
SleepScoreLFP = PickSWTHChannel(basePath,...
                            scoretime,SWWeightsName,...
                            Notch60Hz,NotchUnder3Hz,NotchHVS,NotchTheta,...
                            SWChannels,ThetaChannels,rejectChannels,...
                            overwrite,'noPrompts',noPrompts);

%% CLUSTER STATES BASED ON SLOW WAVE, THETA, EMG

%Calculate the scoring metrics: broadbandLFP, theta, EMG in 
disp('Quantifying metrics for state scoring')
[SleepScoreMetrics,StatePlotMaterials] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,ACCM,overwrite,...
                                           'onSticky',stickytrigger);
                                       
%Use the calculated scoring metrics to divide time into states
disp('Clustering States Based on EMG, SW, and TH LFP channels')
[ints,idx,MinTimeWindowParms] = ClusterStates_DetermineStates(...
                                           SleepScoreMetrics,MinTimeWindowParms);


                                
%% RECORD PARAMETERS from scoring
detectionparms.userinputs = p.Results;
detectionparms.MinTimeWindowParms = MinTimeWindowParms;
detectionparms.SleepScoreMetrics = SleepScoreMetrics;

% note and keep special version of original hists and threshs

SleepState.ints = ints;
SleepState.idx = idx;


SleepState.detectorinfo.detectorname = 'SleepScoreMaster';
SleepState.detectorinfo.detectionparms = detectionparms;
SleepState.detectorinfo.detectionparms.histsandthreshs_orig = detectionparms.SleepScoreMetrics.histsandthreshs;
SleepState.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
SleepState.detectorinfo.StatePlotMaterials = StatePlotMaterials;

%Saving SleepStates
save(bz_sleepstatepath,'SleepState');

%% MAKE THE STATE SCORE OUTPUT FIGURE
%ClusterStates_MakeFigure(stateintervals,stateIDX,figloc,SleepScoreMetrics,StatePlotMaterials);
try
    ClusterStates_MakeFigure(SleepState,basePath,noPrompts);
    disp('Figures Saved to StateScoreFigures')
catch
    disp('Figure making error')
end

%% JOIN STATES INTO EPISODES

% Extract states, Episodes, properly organize params etc, prep for final saving
disp('Calculating/Saving Episodes')
StatesToEpisodes(SleepState,basePath);

disp(['Sleep Score ',recordingname,': Complete!']);

%% PROMPT USER TO MANUALLY CHECK DETECTION WITH THESTATEEDITOR
% if ~noPrompts
%     str = input('Would you like to check detection with TheStateEditor? [Y/N] ','s');
%     switch str
%         case {'Y','y',''}
%             TheStateEditor([basePath,filesep,recordingname])
%         case {'N','n'}
%         otherwise
%             disp('Unknown input..... you''ll have to load TheStateEditor on your own')
%     end
% end

end
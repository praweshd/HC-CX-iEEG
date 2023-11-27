function [SleepScoreMetrics,StatePlotMaterials] = ClusterStates_GetMetrics(...
    basePath,SleepScoreLFP,EMG,ACCM,overwrite,varargin)
%StateID(LFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,WSEpisodes)
 
%% Params
p = inputParser;
addParameter(p,'onSticky',false)
parse(p,varargin{:})
onSticky = p.Results.onSticky; 

%This is the sticky trigger passed through to DetermineStates via histsandthreshs
if onSticky
    stickySW = true; stickyTH=false; stickyEMG=true;
else
    stickySW = false; stickyTH=false; stickyEMG=false;
end
%% Buzcode name of the SleepScoreMetrics.LFP.mat file
[datasetfolder,recordingname,extension] = fileparts(basePath);
recordingname = [recordingname,extension]; % fileparts parses '.' into extension
matfilename = fullfile(basePath,[recordingname,'.SleepScoreMetrics.LFP.mat']);
plotmaterialsfilename = fullfile(basePath,[recordingname,'.StatePlotMaterials.mat']);

if exist(matfilename) & exist(plotmaterialsfilename) & overwrite == false
    load(matfilename,'SleepScoreMetrics')
    load(plotmaterialsfilename,'StatePlotMaterials')
    return
end


%Get the SW weights from SleepScoreLFP - if it's not there, load the
%weights
try
    SWweights = SleepScoreLFP.params.SWweights;
    SWfreqlist = SleepScoreLFP.params.SWfreqlist;
catch
    load('SWweights.mat')
end
%% Downsample and filter the LFP from PickSWTHChannel
%Make Downsample to niquest frequency

if SleepScoreLFP.sf == 1250
    downsamplefactor = 5;
elseif SleepScoreLFP.sf == 250
    downsamplefactor = 1;
elseif SleepScoreLFP.sf == 1000
    downsamplefactor = 4;
else
    display('sf not recognized... if only you made this able to set its own downsample...')
end
swLFP = downsample(SleepScoreLFP.swLFP,downsamplefactor);
thLFP = downsample(SleepScoreLFP.thLFP,downsamplefactor);
sf_LFP = SleepScoreLFP.sf/downsamplefactor;


%% Calculate broadbandslowwave metric
%display('FFT Spectrum for Broadband LFP')
%Timing Parameters
window = 10;   %s
noverlap = 9;  %s
smoothfact = 10; %units of seconds - smoothing factor

if strcmp(SWweights,'PSS')
    %Put the LFP in the right structure format
    lfp.data = swLFP;
    lfp.timestamps = SleepScoreLFP.t;
    lfp.samplingRate = SleepScoreLFP.sf;
    %Calculate PSS
    [specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,'frange',[4 90]);
    broadbandSlowWave = -specslope.data; %So NREM is higher as opposed to lower
    t_clus = specslope.timestamps;
    swFFTfreqs = specslope.freqs;
    swFFTspec = 10.^spec.amp; %To reverse log10 in bz_PowerSpectrumSlope
    badtimes = false;
    %ADD HERE: Bad times detection using swFFTspec similar to below. make bad times nan
   % SWfreqlist = specslope.freqs;
else
    freqlist = logspace(0,2,100);
    [swFFTspec,swFFTfreqs,t_clus] = spectrogram(single(swLFP),window*sf_LFP,noverlap*sf_LFP,freqlist,sf_LFP);
    swFFTspec = abs(swFFTspec);
    [zFFTspec,mu,sig] = zscore(log10(swFFTspec)');
    % Remove transients before calculating SW histogram
    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;

    %Calculate per-bin weights onto SlowWave
    assert(isequal(freqlist,SWfreqlist),...
        'spectrogram freqs.  are not what they should be...')
    broadbandSlowWave = zFFTspec*SWweights';
    
end

%Smooth and 0-1 normalize
broadbandSlowWave = smooth(broadbandSlowWave,smoothfact./mean(diff(t_clus)));
broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));

 
%% Smooth and 0-1 normali
thsmoothfact = 10; %used to be 15

%% Calculate theta
%display('FFT Spectrum for Theta')

% %NarrowbandTheta
f_all = [2 20];
f_theta = [5 10];
freqlist = logspace(log10(f_all(1)),log10(f_all(2)),100);


[thFFTspec,thFFTfreqs,t_thclu] = spectrogram(single(thLFP),window*sf_LFP,noverlap*sf_LFP,freqlist,sf_LFP);
thFFTspec = (abs(thFFTspec));
[~,mu_th,sig_th] = zscore(log10(thFFTspec)');

thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
allpower = sum((thFFTspec),1);
thpower = sum((thFFTspec(thfreqs,:)),1);

thratio = thpower./allpower;    %Narrowband Theta
thratio = smooth(thratio,thsmoothfact./mean(diff(t_thclu)));
thratio = (thratio-min(thratio))./max(thratio-min(thratio));
 
%% EMG
dtEMG = 1/EMG.samplingFrequency;
EMG.smoothed = smooth(EMG.data,smoothfact/dtEMG,'moving');

reclength = round(EMG.timestamps(end)); %What does this get used for?

%interpolate to FFT time points;
%t_EMG = interp1(EMG.timestamps,EMG.timestamps,t_clus,'nearest');% use t_clus
EMG = interp1(EMG.timestamps,EMG.smoothed,t_clus,'nearest');

%Min/Max Normalize
EMG = (EMG-min(EMG))./max(EMG-min(EMG));


%% Divide PC1 for SWS
%Note: can replace all of this with calls to bz_BimodalThresh
numpeaks = 1;
numbins = 12;
%numbins = 12; %for Poster...
while numpeaks ~=2
    [swhist,swhistbins]= hist(broadbandSlowWave,numbins);
    
    [PKS,LOCS] = findpeaks_SleepScore(swhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end


betweenpeaks = swhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-swhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

swthresh = betweenpeaks(diploc);

%Set transients to wake state
broadbandSlowWave(badtimes,1)=swhistbins(LOCS(1));
ACCM = ACCM(1:length(broadbandSlowWave));
 
%PPU INTERVENTION 
%SWS time points
if length(ACCM)~= 0
    disp('NREM: Getting wake times from accelerometer data')
    disp('Size of broadband')
    disp(size(broadbandSlowWave))
    disp('Size of ACCM')
    disp(size(ACCM))
    NREMtimes = (broadbandSlowWave >swthresh & ACCM< 0.5);
elseif length(ACCM) == 0 
    NREMtimes = (broadbandSlowWave >swthresh);
end 


%% Then Divide EMG
numpeaks = 1;
numbins = 12;
if sum(isnan(EMG))>0
   error('EMG seems to have NaN values...') 
end

while numpeaks ~=2
    [EMGhist,EMGhistbins]= hist(EMG(NREMtimes==0),numbins);
    %[EMGhist,EMGhistbins]= hist(EMG,numbins);

    [PKS,LOCS] = findpeaks_SleepScore([0 EMGhist],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    
    if numpeaks ==100
        display('Something is wrong with your EMG')
        return
    end
end

betweenpeaks = EMGhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-EMGhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

EMGthresh = betweenpeaks(diploc);

%PPU INTERVENTION 
if length(ACCM)~= 0
    disp('MOVTIMES: Getting wake times from accelerometer data')
    MOVtimes = (broadbandSlowWave(:)<swthresh & ACCM(:)>0.5);
elseif length(ACCM) == 0 
    MOVtimes = (broadbandSlowWave(:)<swthresh & EMG(:)>EMGthresh);
end 


%% Then Divide Theta
numpeaks = 1;
numbins = 12;
while numpeaks ~=2 && numbins <=25
    %[THhist,THhistbins]= hist(thratio(SWStimes==0 & MOVtimes==0),numbins);
    [THhist,THhistbins]= hist(thratio(MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

numbins = 12;
%numbins = 15; %for Poster...
while numpeaks ~=2 && numbins <=25
    [THhist,THhistbins]= hist(thratio(NREMtimes==0 & MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

if length(PKS)==2
    betweenpeaks = THhistbins(LOCS(1):LOCS(2));
    [dip,diploc] = findpeaks_SleepScore(-THhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

    THthresh = betweenpeaks(diploc);
    %PPU INTERVENTION 
    if length(ACCM)~= 0
        disp('REMTIMES: Getting wake times from accelerometer data')
        REMtimes = (broadbandSlowWave<swthresh & ACCM<0.5 & thratio>THthresh);
    elseif length(ACCM) == 0 
        REMtimes = (broadbandSlowWave<swthresh & EMG<EMGthresh & thratio>THthresh);
    end 
else
    THthresh = 0;
%     REMtimes =(broadbandSlowWave<swthresh & EMG<EMGthresh);
end

histsandthreshs = v2struct(swhist,swhistbins,swthresh,EMGhist,EMGhistbins,...
    EMGthresh,THhist,THhistbins,THthresh,...
    stickySW,stickyTH,stickyEMG);

%% Ouput Structure: StateScoreMetrics
LFPparams = SleepScoreLFP.params;
THchanID = SleepScoreLFP.THchanID; SWchanID = SleepScoreLFP.SWchanID;

SleepScoreMetrics = v2struct(broadbandSlowWave,thratio,EMG,...
    t_clus,badtimes,reclength,histsandthreshs,LFPparams,THchanID,SWchanID,...
    recordingname);
%save(matfilename,'SleepScoreMetrics');

StatePlotMaterials = v2struct(swFFTfreqs,swFFTspec,thFFTfreqs,thFFTspec);
%save(plotmaterialsfilename,'StatePlotMaterials'); 

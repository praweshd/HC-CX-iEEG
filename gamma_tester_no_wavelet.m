%% USAGE - Tester step-by-step 
% For Gamma event detection 
% Prawesh Dahal
% Revised: June 20, 2019                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       , 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================
%=============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Parameters Run This 
close all
clear all
clc

state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5}; 
NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM;
 
num_CH = 205;
gam_CH = 10;   

%Parameters
Rs = 1250; 
resample_factor = 1; 
R_s =Rs/resample_factor;
freq = 1:120; 

wini = R_s; 
winf = 1*R_s; 
smooth_win= round(R_s*(0.01));

arburg_n=2;

lowThresholdFactor = -50; % Ripple envoloppe must exceed lowThresholdFactor*stdev
highThresholdFactor = 0; %

gamb = [30 60];
lowb = [1 16] ;
highb = [80 120]; 

maxg = 300;
ming = 50;
interg = 100; 

tic
lfp_file = dir('*.lfp');
CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));


%% ================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This tests some of the spindle detect tests

%% Obtain States
state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5};
NREM = or(StateIntervals{2}, StateIntervals{3});
% NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM;

% Define thresholds
lowThresholdFactor = 2;
highThresholdFactor = 4;

% Load file and restrict to state
num_CH = 205;
spi_CH = 147;

filename = dir('*.lfp');
[~, fbasename, ~] = fileparts(filename.name);
lfp = LoadLfp(fbasename,num_CH,spi_CH); 
restricted_lfp = [Range(Restrict(lfp,state),'s') Data(Restrict(lfp,state))];

%% %% Downsample 
close all
Fs = 1/125; 

dLfpr = restricted_lfp(:,2);
dLfp = resample(dLfpr, 125, 1250); 
rg = restricted_lfp(:,1);
rg = rg(1:10:end);

time = 0:Fs:(length(dLfp))*Fs;
time = time(2:end); 

figure (1); plot(time, dLfp); title('Step 1: Downsampled to 1250 Hz'); xlabel('time(s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfft = 1024;
fs = 125; 
window = nfft; 
noverlap = nfft/2;

[s,f,t] = spectrogram(dLfp,window,noverlap,[],fs);


figure
spectrogram(dLfp,window,noverlap,[],fs,'yaxis')
colormap jet

%% Whiten
figure (2)
 A = arburg(dLfp,2);
Wdlfp = Filter0(A, dLfp);
plot(time, Wdlfp) 

%% Spectrogram  
 
[s,f,t] = spectrogram(Wdlfp,window,noverlap,[],fs);

figure
spectrogram(Wdlfp,window,noverlap,[],fs,'yaxis')
colormap jet
%% % Bandpass Filter
Fs = 1/125; 
Wn = [25 50]; % spindle freq range [Hz]
rs = 125; % new sampling rate
[b,a] = butter(9, 2*Wn/rs, 'bandpass');
fil_sleep = filtfilt(b,a,dLfp);
fil_sleep = [rg fil_sleep];

% Plot BP signal
Rs = 125;
time = 0:Fs:(length(fil_sleep(:,2)))*Fs;
time = time(2:end); 
FIL =  fil_sleep(:,2);
zoomup_i = 0.5*10^6; 
zoomup_f = zoomup_i + 10*Rs; 
 
figure (2) 
subplot(2,1,1); plot(time, FIL); title(['Step 1: BP ', num2str(Wn(1)),'-', num2str(Wn(2)), ' Hz NREM signal']); xlabel('time (s)');
subplot(2,1,2); plot(time(zoomup_i:zoomup_f), FIL(zoomup_i:zoomup_f)); title('Zoom'); xlabel('time (s)'); 

[s,f,t] = spectrogram(FIL,window,noverlap,[],fs);

figure
spectrogram(FIL,window,noverlap,[],fs,'yaxis')
colormap jet
%% Square
  
FILSQ = FIL.^2;
figure (3) 
plot(time(zoomup_i:zoomup_f), FILSQ(zoomup_i:zoomup_f)); title('Step 2: Squared'); xlabel('time (s)'); 

%% Window
windowLength = round(125/1250*11);
window = ones(windowLength,1)/windowLength;

FILWIN = Filter0(window,sum(FILSQ,2)); 

figure (4) 
plot(time(zoomup_i:zoomup_f), FILWIN(zoomup_i:zoomup_f))
hold on
plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*mean(FILWIN),'r')
title('Step 3: Window'); xlabel('time (s)'); 

%% Normalize  
sd = []; keep = []; restrict =[];  
keep = [];
if ~isempty(restrict),
	keep = fil_sleep(:,1)>=restrict(1)&fil_sleep(:,1)<=restrict(2);
end

%% Get standard deviation and normalize 
[~,sd] = unity(FILWIN,sd,keep);
disp(sd)

[normalizedSquaredSignal,sd] = unity(FILWIN,sd,keep);

figure (5) 
subplot(2,1,1)
plot(time, normalizedSquaredSignal)
hold on
plot(time, ones(length(time),1)*mean(normalizedSquaredSignal), 'r') 

subplot(2,1,2)
plot(time(zoomup_i:zoomup_f), normalizedSquaredSignal(zoomup_i:zoomup_f))
hold on 
plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*mean(normalizedSquaredSignal),'r')
title('Step 4: Normalized signal - mean = 0, std = 1'); 
xlabel('time (s)'); xlim([4000 4010]); %ylim([-0.25 2]);

%% Threshold (FIRST-PASS)
 
thresholded = normalizedSquaredSignal > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);

% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1,
	start = start(1:end-1);
end
% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start),
    stop = stop(2:end);
end
% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1),
	stop(1) = [];
	start(end) = [];
end

firstPass = [start,stop];

%% SECOND PASS
minInterRippleInterval = 30; 

% Merge ripples if inter-ripple period is too short
minInterRippleSamples = minInterRippleInterval/1000*Rs;
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - ripple(2) < minInterRippleSamples,
		% Merge
		ripple = [ripple(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; ripple];
		ripple = firstPass(i,:);
	end
end
secondPass = [secondPass ; ripple];

%% THIRD PASS

% Discard ripples with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
	[maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
	if maxValue > highThresholdFactor,
		thirdPass = [thirdPass ; secondPass(i,:)];
		peakNormalizedPower = [peakNormalizedPower ; maxValue];
	end
end

%%
 
fig=figure_ctrl('SPI Detect',1500,1000);
subplot(3,1,1)
plot(time(zoomup_i:zoomup_f), normalizedSquaredSignal(zoomup_i:zoomup_f))
hold on 
plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*mean(normalizedSquaredSignal),'r')
title('Normalized signal - mean = 0, std = 1'); 
xlabel('time (s)'); 
ylim([-0.25 3.5]);

subplot(3,1,2)
% plot(time(1:zoomup), thresholded(1:zoomup))
plot(time(zoomup_i:zoomup_f),thresholded(zoomup_i:zoomup_f))

subplot(3,1,3)
% plot(time(2:zoomup), diff(thresholded(1:zoomup)))
plot(time(zoomup_i+1:zoomup_f),diff(thresholded(zoomup_i:zoomup_f)))
 
%% Final Flow Plot 
fig=figure_ctrl('SPI Detect Step',1500,1000);
picnum = 4; 
subplot(picnum,1,1)
plot(time(zoomup_i+1:zoomup_f), dLfp(zoomup_i+1:zoomup_f))

subplot(picnum,1,2)
plot(time(zoomup_i+1:zoomup_f), FIL(zoomup_i+1:zoomup_f))

subplot(picnum,1,3)
plot(time(zoomup_i+1:zoomup_f), FILSQ(zoomup_i+1:zoomup_f))

subplot(picnum,1,4)
plot(time(zoomup_i+1:zoomup_f), normalizedSquaredSignal(zoomup_i+1:zoomup_f))
hold on 
plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*mean(normalizedSquaredSignal),'r')
plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*lowThresholdFactor,'g')
plot(time(zoomup_i:zoomup_f), ones(length(time(zoomup_i:zoomup_f)),1)*highThresholdFactor,'c')
xlim([4000 4010]); ylim([-0.5 8]);

%%
fig=figure_ctrl('SPI Firstpass',1500,1000);
plot(time(zoomup_i+1:zoomup_f), normalizedSquaredSignal(zoomup_i+1:zoomup_f))
hold on 
plot(startT, ones(length(startT))*2,'*')
% plot(stopT, ones(length(stopT))*2,'*')
xlim([4000 4010]); ylim([-0.5 8]);

% subplot(picnum,1,5)
% plot(time(zoomup_i+1:zoomup_f), dLfp(zoomup_i+1:zoomup_f))
% 
% subplot(picnum,1,6)
% plot(time(zoomup_i+1:zoomup_f), dLfp(zoomup_i+1:zoomup_f))
% 
 
%%

function y = Filter0(b,x)

if size(x,1) == 1,
	x = x(:);
end

if mod(length(b),2) ~= 1,
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

end 

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;

end 




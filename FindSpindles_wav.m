function [spindles] = FindSpindles_wav(filtered,varargin)
%USAGE 
% To detect spindle events from wavelet ratio file
% Modified - Prawesh Dahal/Naureen Ghani (BRAIN paper)
% Inspired from: 
% FindRipples - Find hippocampal Ripples (100~200Hz oscillations).
%
%  USAGE
%
%    [spindles,stdev,noise] = Findspindles(filtered,<options>)
%
%    spindles are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.
%
%    filtered       spindle-band filtered LFP (one channel).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for spindle beginning/end and peak, in multiples
%                   of the stdev (default = [2 5])
%     'durations'   min inter-spindle interval, max spindle duration and min spindle duration, in ms
%                   (default = [30 100 20])
%     'baseline'    interval used to compute normalization (default = all)
%     'restrict'    same as 'baseline' (for backwards compatibility)
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy spindle-band filtered channel used to exclude spindle-
%                   like noise (events also present on this channel are
%                   discarded)
%    =========================================================================
%
%  OUTPUT
%
%    spindles        for each spindle, [start_t peak_t end_t peakNormalizedPower]
%    stdev          standard deviation of the NSS (can be reused subsequently)
%    noise          spindle-like activity recorded simultaneously on the noise
%                   channel (for debugging info)
%
%  SEE
%
%    See also FilterLFP, spindlestats, SavespindleEvents, Plotspindlestats.

% Copyright (C) 2004-2011 by Michaël Zugaro, initial algorithm by Hajime Hirase
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
frequency = 125;
show = 'off';
restrict = [];
sd = [];
lowThresholdFactor = -0.1; % spindle envoloppe must exceed lowThresholdFactor*stdev
highThresholdFactor = 0; % spindle peak must exceed highThresholdFactor*stdev
minInterspindleInterval = 400; % in ms
maxspindleDuration = 3000; % in ms
minspindleDuration = 300; % in ms
noise = [];
pow_thresh = 150;

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(filtered) | size(filtered,2) ~= 2,
	error('Parameter ''filtered'' is not a Nx2 matrix (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
% 	if ~ischar(varargin{i}),
% 		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).']);
% 	end
	switch(lower(varargin{i})),
		case 'thresholds',
			thresholds = varargin{i+1};
% 			if ~isivector(thresholds,'#2','>0'),
% 				error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
% 			end
			lowThresholdFactor = thresholds(1);
			highThresholdFactor = thresholds(2);
		case 'durations',
			durations = varargin{i+1};
% 			if ~isivector(durations,'#3','>0'),
% 				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
% 			end
			minInterspindleInterval = durations(1);
			maxspindleDuration = durations(2);
			minspindleDuration = durations(3);
        case 'frequency',
			frequency = varargin{i+1};
% 			if ~isdscalar(frequency,'>0'),
% 				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
% 			end
		case 'show',
			show = varargin{i+1};
% 			if ~isstring(show,'on','off'),
% 				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
% 			end
		case {'baseline','restrict'},
			restrict = varargin{i+1};
% 			if ~isempty(restrict) & ~isdvector(restrict,'#2','<'),
% 				error('Incorrect value for property ''restrict'' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
% 			end
		case 'stdev',
			sd = varargin{i+1};
% 			if ~isdscalar(sd,'>0'),
% 				error('Incorrect value for property ''stdev'' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
% 			end
		case 'noise',
			noise = varargin{i+1};
% 			if ~isdmatrix(noise) | size(noise,1) ~= size(filtered,1) | size(noise,2) ~= 2,
% 				error('Incorrect value for property ''noise'' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).');
% 			end
        case 'power_thresh',
			pow_thresh = varargin{i+1};
% 		otherwise,
% 			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Findspindles">Findspindles</a>'' for details).']);
	end
end

% % Parameters
% windowLength = round(frequency/1250*11);
% 
% % Square and normalize signal
% signal = filtered(:,2);
% squaredSignal = signal.^2;
% window = ones(windowLength,1)/windowLength;
% keep = [];
% if ~isempty(restrict),
% 	keep = filtered(:,1)>=restrict(1)&filtered(:,1)<=restrict(2);
% end
% 
% [normalizedSquaredSignal,sd] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);

normalizedSquaredSignal = filtered(:, 2);

% Detect spindle periods by thresholding normalized squared signal
thresholded = normalizedSquaredSignal > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last spindle if it is incomplete
if isempty(start),
	disp('Detection by thresholding failed');
    spindles = [0 0 0 0];
	return
else
    if length(stop) == length(start)-1,
        start = start(1:end-1);
    end
    % Exclude first spindle if it is incomplete
    if length(stop)-1 == length(start),
        stop = stop(2:end);
    end
    % Correct special case when both first and last spindles are incomplete
    if start(1) > stop(1),
        stop(1) = [];
        start(end) = [];
    end
    firstPass = [start,stop];
    disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

% Merge spindles if inter-spindle period is too short
minInterspindlesamples = minInterspindleInterval/1000*frequency;
secondPass = [];
spindle = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - spindle(2) < minInterspindlesamples,
		% Merge
		spindle = [spindle(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; spindle];
		spindle = firstPass(i,:);
	end
end
secondPass = [secondPass ; spindle];
if isempty(secondPass),
	disp('spindle merge failed');
    spindles = [0 0 0 0];
	return
else
	disp(['After spindle merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard spindles with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
	[maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
	if maxValue > highThresholdFactor,
		thirdPass = [thirdPass ; secondPass(i,:)];
		peakNormalizedPower = [peakNormalizedPower ; maxValue];
	end
end
if isempty(thirdPass),
	disp('Peak thresholding failed.');
    spindles = [0 0 0 0];
	return
else
	disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end
% 
% % Detect negative peak position for each spindle
% peakPosition = zeros(size(thirdPass,1),1);
% for i=1:size(thirdPass,1),
% 	[minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
% 	peakPosition(i) = minIndex + thirdPass(i,1) - 1;
% end

% Detect peak position for each spindle
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1),
	[minValue,minIndex] = max(normalizedSquaredSignal(thirdPass(i,1):thirdPass(i,2)));
	peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end

% Discard spindles that are way too long
time = filtered(:,1);
spindles = [time(thirdPass(:,1)) time(peakPosition) time(thirdPass(:,2)) peakNormalizedPower];
duration = spindles(:,3)-spindles(:,1);
spindles(duration>maxspindleDuration/1000,:) = [];
disp(['After max duration test: ' num2str(size(spindles,1)) ' events.']);

%Discard spindles that are way too short
duration = spindles(:,3)-spindles(:,1);
spindles(duration<minspindleDuration/1000,:) = [];
disp(['After min duration test: ' num2str(size(spindles,1)) ' events.']);

% Exclude spindle if it bridges a movement interval
fourthPass = spindles;
mvmt_breaks = find(diff(filtered(:, 1))>(1.9/frequency));
time_breaks = filtered(mvmt_breaks, 1);
for ii = 1:size(fourthPass, 1)
    for jj = 1:size(time_breaks, 1)
       if fourthPass(ii, 1) <= time_breaks(jj) && fourthPass(ii, 3) >= time_breaks(jj)
           spindles(ii, :) = [0 0 0 0];
       end
    end
end
spindles_only_index = find(spindles(:, 1) ~= 0);
spindles = spindles(spindles_only_index, :);

% Exclude spindle if it has too low power
low_index = find(spindles(:, 4)<pow_thresh);
spindles(low_index, :) = [];

% If a noisy channel was provided, find spindle-like events and exclude them
% bad = [];
% if ~isempty(noise),
% 	% Square and pseudo-normalize (divide by signal stdev) noise
% 	squaredNoise = noise(:,2).^2;
% 	window = ones(windowLength,1)/windowLength;
% 	normalizedSquaredNoise = unity(Filter0(window,sum(squaredNoise,2)),sd,[]);
% 	excluded = logical(zeros(size(spindles,1),1));
% 	% Exclude spindles when concomittent noise crosses high detection threshold
% 	previous = 1;
% 	for i = 1:size(spindles,1),
% 		j = FindInInterval(noise,[spindles(i,1),spindles(i,3)],previous);
% 		previous = j(2);
% 		if any(normalizedSquaredNoise(j(1):j(2))>(highThresholdFactor)),            %Changed to highThresholdFactor-2 to increase sensitivity of noise channel
% 			excluded(i) = 1;
% 		end
% 	end
% 	bad = spindles(excluded,:);
% 	spindles = spindles(~excluded,:);
% 	disp(['After noise removal: ' num2str(size(spindles,1)) ' events.']);
% end

% Optionally, plot results
if strcmp(show,'on'),
	figure;
	if ~isempty(noise),
		MultiPlotXY([time signal],[time squaredSignal],[time normalizedSquaredSignal],[time noise(:,2)],[time squaredNoise],[time normalizedSquaredNoise]);
		nPlots = 6;
		subplot(nPlots,1,3);
 		ylim([0 highThresholdFactor*1.1]);
		subplot(nPlots,1,6);
  		ylim([0 highThresholdFactor*1.1]);
	else
		MultiPlotXY([time signal],[time squaredSignal],[time normalizedSquaredSignal]);
%  		MultiPlotXY(time,signal,time,squaredSignal,time,normalizedSquaredSignal);
		nPlots = 3;
		subplot(nPlots,1,3);
  		ylim([0 highThresholdFactor*1.1]);
	end
	for i = 1:nPlots,
		subplot(nPlots,1,i);
		hold on;
  		yLim = ylim;
		for j=1:size(spindles,1),
			plot([spindles(j,1) spindles(j,1)],yLim,'g-');
			plot([spindles(j,2) spindles(j,2)],yLim,'k-');
			plot([spindles(j,3) spindles(j,3)],yLim,'r-');
			if i == 3,
				plot([spindles(j,1) spindles(j,3)],[spindles(j,4) spindles(j,4)],'k-');
			end
		end
		for j=1:size(bad,1),
			plot([bad(j,1) bad(j,1)],yLim,'k-');
			plot([bad(j,2) bad(j,2)],yLim,'k-');
			plot([bad(j,3) bad(j,3)],yLim,'k-');
			if i == 3,
				plot([bad(j,1) bad(j,3)],[bad(j,4) bad(j,4)],'k-');
			end
		end
		if mod(i,3) == 0,
			plot(xlim,[lowThresholdFactor lowThresholdFactor],'k','linestyle','--');
			plot(xlim,[highThresholdFactor highThresholdFactor],'k-');
		end
	end
end

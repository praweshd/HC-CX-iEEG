function [ ripples, sd ] = RIP_Detect_Cx( num_CH, rip_CH, noise_CH, varargin )
%USAGE 
%Ripple Detection
% 
% input: total number of channels in recording, channel for ripple analysis
%          (+1 from neuroscope channel number), ensure .lfp file and -states.mat(from StateEditor) is in the folder
%      : uses OOP to load lfp and restrict to states
%      : uses FMAtoolbox adapted for detection, file saving, event file
%        saving
% ripple band = [100 250]
% durations = [30 100 20] (min inter-ripple interval, max ripple duration and
%                           min ripple duration, in ms)
% output: .mat file (named with fbasename+rip_CH)
%          column 1 = start; column 2 = peak; column 3 = end; column 4 =
%          power
%         .rip.evt file 


%=========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'rip_thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5])
%     'state'           select state between NREM, REM, WAKE
%     'baseline'        select baseline for calculation of standard
%                       deviation (as a vector of start and stop in
%                       seconds)
%     'stdev'           std deviation to use for calculation

%Defaults


% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'rip_thresholds',
			rip_thresholds = varargin{i+1};
			if ~isivector(rip_thresholds,'#2','>0'),
				error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
			end
			rip_lowThresholdFactor = rip_thresholds(1);
			rip_highThresholdFactor = rip_thresholds(2);
        case 'state',
			state = varargin{i+1};
            state_name = state;
			if ~isstring(state,'NREM','REM', 'WAKE'),
				error('Incorrect value for property ''state'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'baseline',
			baseline = varargin{i+1};
			if ~isivector(baseline,'#2','>0'),
				error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
            end
        case 'stdev',
			sd = varargin{i+1};
			if ~isdscalar(sd,'>0'),
				error('Incorrect value for property ''stdev'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
			end
    end
end
            
filename = dir('*.lfp');
[pathstr, fbasename, fileSuffix] = fileparts(filename.name);
nchannels = num_CH;
ripple_channel = rip_CH;  
noise_channel = noise_CH;
channelID = ripple_channel - 1;

if ~isempty(dir('*-newstates*'))
    state_mat = dir('*-newstates*');
    disp('Loading CLEANED states file')
    load (state_mat.name);
    states = newstates;
else
    state_mat = dir('*-states*');
    load (state_mat.name);
end
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5};
NREM = StateIntervals{3};
WAKE = StateIntervals{1};

% State parameter
state = NREM;
% if strcmp(state,'NREM'),
%     state = NREM;
% elseif strcmp(state, 'REM'),
%         state = REM;
% else strcmp(state, 'WAKE'),
%     state = WAKE;
% end

Rs=1250;
disp('Loading lfp')
lfp = LoadLfp(fbasename,nchannels,ripple_channel);  
sleep=[Range(Restrict(lfp, state), 's') Data(Restrict(lfp, state))];

signal=sleep(:,2);
n=3;
Wn=[110 180]; %[110 180]; %% [110 180]; %% [85 150]
[b,a]=butter(n,2*Wn/Rs,'bandpass'); 

fil_sleep_V=filtfilt (b,a,signal);
fil_sleep=[sleep(:,1) fil_sleep_V];

%fil_sleep = FilterLFP([Range(Restrict(lfp, state), 's') Data(Restrict(lfp, state))], 'passband', [100 250]);

if ~isempty(noise_CH)
    disp(['Using noise CH ', num2str(noise_CH)])
    noise_lfp = LoadLfp(fbasename,nchannels,noise_channel); 
    noise= [Range(Restrict(noise_lfp, state), 's') Data(Restrict(noise_lfp, state))];
    fil_noise_V=filtfilt (b, a,noise(:,2));
    fil_noise=[noise(:,1) fil_noise_V];    
else 
    disp('No noise CH')
    fil_noise = [];    
end

% fil_noise = FilterLFP([Range(Restrict(noise_lfp, state), 's') Data(Restrict(noise_lfp, state))], 'passband', [100 250]);
% [ripples,sd,bad] = FindRipples(fil_sleep,'thresholds', [rip_lowThresholdFactor rip_highThresholdFactor], 'durations', [30 110 30], 'noise', fil_noise); % 'stdev', sd, 'baseline', baseline, 
[ripples,sd,bad] = FindRipples(fil_sleep,'thresholds', [rip_lowThresholdFactor rip_highThresholdFactor], 'durations', [30 90 20], 'noise', fil_noise); % 'stdev', sd, 'baseline', baseline, 

% ripple_file = strcat(fbasename,'_CH', num2str(rip_CH),'24_85_150', state_name);
% save (ripple_file, 'ripples')
% ripple_events = strcat(fbasename, num2str(rip_CH), state_name, '24_85_150.cxd.evt');
% SaveRippleEvents(ripple_events,ripples,channelID);

end
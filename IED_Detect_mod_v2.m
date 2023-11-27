function IED = IED_Detect_mod_v2(num_CH, IED_CH, varargin)
% USAGE: IED_DETECT Detects IED events in LFP file
%       IED = IED_DETECT(NUM_CH, IED_CH) uses default parameters
%       IED = IED_DETECT(NUM_CH, IED_CH, passband, state, med_thresh,
%       ied_int, minburstInt, mean_thresh) uses specified parameters 
%
% REQUIREMENTS
%       NUM_CH: total number of channels
%       IED_CH: channel for IED analysis
%       Note:   must be in path of .lfp file and -states.mat(from StateEditor)
%                   file
%
% OPTIONAL INPUT ARGUMENTS
%       PASSBAND:    specifies a range for passband filtering 
%       STATE:       select between WAKE, NREM, or REM
%       MED_THRESH:  specifies a median threshold for enhancing LFP filter
%       IED_INT:     specifies a time interval between IED events 
%       MINBURSTINT: specifies a time interval between IED bursts
%       MEAN_THRESH: specifies a mean threshold for unfiltered LFP
%
% OUTPUT
%       IED.mat
%           named with fbasename+IED_CH where column 1 = start; 
%           column 2 = peak; column 3 = end; column 4 = power;
%       .ied.evt file
%           named with fbasename+IED_CH
%
% EXAMPLES:
% [IED] = IED_Detect(128,94)
% [IED] = IED_DETECT(128,94, [25 80], 'NREM', 5, 0.2, 1, 2);
%
% Prawesh Dahal and Jennifer Gelinas (2018) 

%Parse Inputs
p = inputParser;
defaultPassband = [25 80];
defaultState = 'NREM';
defaultMed_thresh = 5;
defaultIed_int = 0.2;
defaultMinburstInt = 0.1;
defaultMeanthresh = 2.25;

addRequired(p, 'num_CH', @isnumeric);
addRequired(p, 'IED_CH', @isnumeric);
addOptional(p, 'passband', defaultPassband, @isivector);
addOptional(p, 'state', defaultState, @is_string);
addOptional(p, 'med_thresh', defaultMed_thresh, @isdscalar);
addOptional(p, 'ied_int', defaultIed_int, @isdscalar);
addOptional(p, 'minburstInt', defaultMinburstInt, @isdscalar);
addOptional(p, 'mean_thresh', defaultMeanthresh, @isdscalar);

parse(p, num_CH, IED_CH , varargin{:});

% Obtain States
state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5};
NREM = or(StateIntervals{2}, StateIntervals{3});
WAKE = StateIntervals{1};

% Define state
state = p.Results.state;

% State parameter
if strcmp(state,'NREM')
    state = NREM;
elseif strcmp(state, 'REM')
    state = REM;
else
    strcmp(state,'WAKE')
    state = WAKE;
end

% Obtain File Name
filename = dir('*.lfp');
[~, fbasename, ~] = fileparts(filename.name);

% Define Passband
IED_lowpassband = p.Results.passband(1);
IED_highpassband = p.Results.passband(2);

% Filter LFP as per passband range provided
lfp = LoadLfp(fbasename,num_CH,IED_CH);
fil_sleep = FilterLFP([Range(Restrict(lfp, state), 's') Data(Restrict(lfp, state))], 'passband', [IED_lowpassband IED_highpassband]);

% Obtain NREM duration to use as a parameter in peak detection
[~,ind] = find(states==3);
transition = cell(1);
for i = 1:length(ind)-1
    if ind(i+1)-ind(i) ~= 1
        transition{i} = [ind(i) ind(i+1)];
    else      
    end
end
transition = transition(~cellfun(@isempty, transition));
transition = horzcat(ind(1),cell2mat(transition),ind(end));

NREM_start = transition(1:2:end);
NREM_end = transition(2:2:end);
[~,ind] = max(NREM_end-NREM_start);
start_NREM = NREM_start(ind);
end_NREM = NREM_end(ind); 

% Define parameters
med_thresh = p.Results.med_thresh;
ied_int = p.Results.ied_int;
% minburstInt = p.Results.minburstInt; 
mean_thresh = p.Results.mean_thresh;

% Detect Peaks in filtered LFP relative to median of filtered LFP
rect_filsleep = abs(fil_sleep);
local_peaks = autofindpeaks(rect_filsleep(:, 1), rect_filsleep(:, 2), 0, median(rect_filsleep(start_NREM:end_NREM, 2))*med_thresh, 100, 100, 1); 

% Sort Local Peaks by Position
[~, I] = sort(local_peaks(:,2)); 
local_peaks = local_peaks(I,:);
% local_peaks is an nx5 matrix such that:
%       Col 1: Position in y-axis
%       Col 2: Position in x-axis
%       Col 3: Height of peak
%       Col 4: Width of peak
%       Col 5: Area of peak
IED_sep = local_peaks;


%%%% Eliminate local peaks based on neurobiology %%%%%

%Step 1: Remove peaks that are other oscillatory events (i.e. ripples)
%based on filtered LFP amplitude (>15Hz)
Rs=1250;
res = floor(IED_sep(:, 2)*Rs);
neg_elim = find(res < 0);
res(neg_elim) = [];

if isempty(res)==1 
       disp('No IED found'); 
       IED = zeros(1,3); 
       IED_tmp = zeros(1,1); 
elseif res == 0
       disp('No IED found'); 
       IED = zeros(1,3); 
       IED_tmp = zeros(1,1); 
elseif length(res) == 1
       disp('No IED found'); 
       IED = zeros(1,3);
       IED_tmp = zeros(1,1); 
else
        lfp2 = Data(lfp);
        Wn = 15;
        [b,a]=butter(3,2*Wn/Rs, 'high');
        mean_lfp = mean(abs(lfp2(start_NREM:end_NREM)));
        frame_time=250e-3; 
        N=round (frame_time * Rs);
        T=res(1:end);

        %data=Dat_tracker(filename.name,T,N,num_CH);
        %data_CHsel = lfp2(:, IED_CH);
        data_CHsel = filtfilt(b, a, lfp2);
         large_IEDs = zeros();
        for ii = 1:length(res)-1
            IED_CHsel(:, ii) = data_CHsel(res(ii)-round(N/2):res(ii)+round(N/2));
            temp = abs(IED_CHsel(:, ii)) > mean_lfp*mean_thresh;
            [~, peak_index]=min (IED_CHsel(:, ii));
            IED_peak_time (ii)= T(ii)-N/2 + peak_index;
            if sum(temp) > 10
                large_IEDs(ii) = 1;
            else
                large_IEDs(ii) = 0;
            end
        end
%         
%         if isempty(data)==1
%             IED = zeros(1,3);
%         else
%         data_CHsel = data(IED_CH, :, :);
%         data_CHsel = squeeze(data_CHsel);
%         data_CHsel = filtfilt(b, a, data_CHsel);
%         large_IEDs = zeros();
%         for ii = 1:length(res)
%             temp = abs(data_CHsel(:, ii)) > mean_lfp*mean_thresh;
%             [~, peak_index]=min (data_CHsel(:, ii));
%             IED_peak_time (ii)= T(ii)-N/2 + peak_index;
%             if sum(temp) > 10
%                 large_IEDs(ii) = 1;
%             else
%                 large_IEDs(ii) = 0;
%             end
%         end

        large_IEDindex = large_IEDs == 1;
        IED_peak_time= IED_peak_time(large_IEDindex);
        IED_tmp=IED_peak_time';

       

end



if length(IED_tmp)>2
        refreactory=0.150*Rs;
        res_diff= diff(IED_tmp);
        res_refract=res_diff>refreactory;
        res_tmp=IED_tmp(2:end).*res_refract;

        %% extraction
        IED= res_tmp(find(res_tmp~=0)); % spikes peak time
    else
        IED=IED_tmp;

end
        IED=IED./Rs;

if isempty(IED)==1 
    disp('No IED Events');
else
    %%% Step 5: Save IED Events File %%%%%
    large_file = strcat(fbasename,'_Ch', num2str(IED_CH), '33_IED');
    save(large_file,'IED');
    IEDs_name = strcat(fbasename, '_', num2str(IED_CH),'_33.ied.evt');
    channelID = IED_CH - 1;

    IEDevts = [IED-0.0625 IED IED+0.0625];
    SaveRippleEvents(IEDs_name, IEDevts, channelID);
    
end




function [ TF_data_crip] = TF_batch_corticalripple( res_CH, num_CH, test_CH, varargin )
% Analysis for cortical ripples (PNAS 2023) 

% Check number of inputs
if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

% Check varargin
if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters  (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
	switch(lower(varargin{i})),
            case 'session'
			session = varargin{i+1};
            case 'analysis'
			analysis = varargin{i+1};
            case 'res_type'
			res_type = varargin{i+1};
            case 'state',
			state = varargin{i+1};
            state_name = state;
			if ~isstring(state,'NREM','REM', 'WAKE'),
				error('Incorrect value for property ''state'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end

if ~isempty(session)
  cd(session);
end

filename = dir('*.lfp');
[pathstr, fbasename, fileSuffix] = fileparts(filename.name);

state_mat = dir('*-states*');
load (state_mat.name);

StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5};
NREM = or(StateIntervals{2}, StateIntervals{3});                            % merges drowsy, NREM, and intermediate sleep
WAKE = StateIntervals{1};
NREM_only = StateIntervals{3};

% State parameter
if strcmp(state,'NREM'),
    state = NREM;
elseif strcmp(state, 'REM'),
        state = REM;
else strcmp(state, 'WAKE'),
    state = WAKE;
end

Rs = 1250;
frame_time = 1;             
offset = 0;             
frequency_vec = 20:1:500; 

xaxis = -0.46:1/Rs:0.46;                                          % 0.175:1/Rs:2.335;   -0.36:1/Rs:1.8
yaxis = 50:1:500;
    
% Calculate the baseline for subsequent z-scoring
k = 1;
for CH = [test_CH]
lfp = LoadLfp(fbasename,num_CH,CH);    
sleep = [Range(Restrict(lfp, state), 's') Data(Restrict(lfp, state))];
    if length(sleep)>500*Rs
        sleep_subset = sleep(1:500*Rs,2);
    else 
        sleep_subset = sleep(:, 2);
    end
    
[S f]=awt_freqlist(sleep_subset,Rs,frequency_vec,'Gabor');              %1:1:500 (ripples); 1:0.5:50 (spindles)

TF_data_crip(k).CH = test_CH(k);
TF_data_crip(k).power_mean = mean(abs(S));
TF_data_crip(k).power_std = std(abs(S));

% Generate 'S' for given res file

if strcmp(res_type, 'ripple') == 1
    res_file = (strcat(fbasename, num2str(res_CH), state_name, 'ripples.mat'));        
end

 [ S,f, wavavg_trials, var_name ] = Wav_res( res_file, num_CH, CH, frame_time, offset, frequency_vec );

% Z-score 'S'

test = TF_data_crip(k).power_mean;
test1 = TF_data_crip(k).power_std;

[a b c] = size(S);

for jj = 1:c
    for ii = 1:b
       S_corr(:, ii, jj) = (abs(S(:, ii, jj))-test(ii))./test1(ii); 
    end
end

% Calculate the average spectrogram and plot the figure

mean_corr = squeeze(mean(S_corr, 3));
figure1 = figure('visible', 'off');
imagesc (xaxis, yaxis, mean_corr(50:1200, 31:481)')          % 800:3500
axis xy;
colormap jet;

print(figure1,'-djpeg ',strcat(fbasename, state_name, num2str(test_CH(k)), '_TF'));
%saveas(figure1, strcat(fbasename, state_name, num2str(test_CH(k)), '_TF', '.fig'));

% Calculate the power spectrum pre/post and plot the figure

S_corr_spec1 = squeeze(mean(S_corr(100:200, :, :), 1));                % 800:2675   500:2375  1250:1875 (1.5s duration)
test2 = mean(S_corr_spec1, 2);
test3 = std(S_corr_spec1, 0, 2);
test4 = test3./sqrt(c);
figure2 = figure('visible', 'off');
shadedErrorBar(yaxis,test2(31:481)',test4(31:481)','g');

hold on

S_corr_spec2 = squeeze(mean(S_corr(625:725, :, :), 1));                % 800:2675   500:2375  1250:1875 (1.5s duration)
test2a = mean(S_corr_spec2, 2);
test3a = std(S_corr_spec2, 0, 2);
test4a = test3a./sqrt(c);
shadedErrorBar(yaxis,test2a(31:481)',test4a(31:481)','b');

print(figure2,'-djpeg ',strcat(fbasename, state_name, num2str(test_CH(k)), '_powspec'));
%saveas(figure2, strcat(fbasename, state_name, num2str(test_CH(k)), '_powspec', '.fig'));

% Calculate the difference in ripple band power before and after

sum_rip_pre = sum(S_corr_spec1(81:231, :));
sum_rip_post = sum(S_corr_spec2(81:231, :));
diff_rip = (sum_rip_post - sum_rip_pre);

TF_data_crip(k).sumrippre = sum_rip_pre;
TF_data_crip(k).sumrippost = sum_rip_post;
TF_data_crip(k).diffspi = diff_rip;

% Calculate the ripple band power over time with ste

sum_rippow = squeeze(sum(S_corr(:, 81:231, :), 2));
mean_rippow = mean(sum_rippow, 2);
ste_rippow = ste(sum_rippow');
TF_data_crip(k).mean_rippow = mean_rippow;
TF_data_crip(k).ste_rippow = ste_rippow;

% Calculate the gamma band power over time with ste

sum_gampow = squeeze(sum(S_corr(:, 41:81, :), 2));
mean_gampow = mean(sum_gampow, 2);
ste_gampow = ste(sum_gampow');
TF_data_crip(k).mean_gampow = mean_gampow;
TF_data_crip(k).ste_gampow = ste_gampow;

k = k+1;
end

% Save the array and copy it

TFdata_file = strcat(fbasename, state_name, res_type, '_TFdata_crip');
save(TFdata_file, 'TF_data_crip');

if ~isempty(analysis)
  copyfile (strcat(TFdata_file, '.mat'), analysis);
end


end


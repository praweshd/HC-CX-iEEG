function [ S,f, wavavg_trials, var_name ] = Wav_res( res_file, num_CH, test_CH, frame_time, offset, frequency_vec )
% Extracts lfp at each timepoint of res_file and calculates wavelet over
% pre-specified time-interval
%   res_file: .mat containing vector of timepoints in first column
%   num_CH: number of channels
%   test_CH: channel selected for analysis
%   frame_time: duration, in s, of interval on which to calculate wavelet;
%   timepoint in res_file will be taken as the centre of frame_time
%   offset: amount of time, in s, to shift centre of frame_time (positive =
%   shift to left, negative = shift to right)
%   frequency_vector: frequencies for use in wavelet calculation (Gabor),
%   in format [low:interval:high]

load(res_file);
var_name = who;
var_name2 = eval(sprintf(var_name{6}));

filename = dir('*.lfp');
[pathstr, fbasename, fileSuffix] = fileparts(filename.name);

% extract the lfp at each timepoint of the res_file
Rs = 1250;
if length(var_name2) > 1000
    res = floor (var_name2(1:1000, 1) *Rs);
else
    res = floor (var_name2(:, 1) *Rs);
end% 1:100 for ripples
                                      
N = frame_time * Rs;
off = offset * Rs;                 
T=res(1:end);                          

data=Dat_tracker(filename.name, T+off, N, num_CH);                                
data_CHsel = squeeze ( data(test_CH, :, :)  );

% calculate the wavelet 
 for ii = 1:length(T)
    [S(:, :, ii),f,psi_array] = awt_freqlist(data_CHsel(:, ii), Rs, frequency_vec,'Gabor');
 end
 
 wavavg_trials = mean(abs(S),3);
 
% wav_file = strcat(fbasename, num2str(test_CH), 'REM', 'wav');
% wavelet.S = S;
% wavelet.f = f;
% %save (wav_file, 'wavelet')
% wavavg_file = strcat(fbasename, num2str(test_CH), var_name{1}, '_wavavg');
% %save (wavavg_file, 'wavavg_trials')
 
%  figure1 = figure;
%  new_xaxis = (-frame_time/2+offset):1/Rs:((frame_time/2 + offset)-1/Rs);
%  new_yaxis = 8:0.5:25;
%  imagesc (new_xaxis, new_yaxis, wavavg_trials(:, 15:49)')
%      axis xy
     %xlim([ 0 N ])
%      hold on
%      line([N/2 N/2], [min(f) max(f)])
     %print(figure1,'-djpeg',cat (2, 'Avgwav',num2str(test_CH), var_name{1}));
end


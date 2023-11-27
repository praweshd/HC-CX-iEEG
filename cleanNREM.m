%% NREM Artifacts Cleaning 
% Prawesh Dahal
% Feb 27, 2020

function cleanNREM(goodch,NREM_start,NREM_end)
%USAGE - NREM artifacts cleaning program
% 
% inputs:
%goodch: one good channel from Neuroscope (matlab-indexed)
%NREM_start: time in seconds (from NS) of start of ~ 5s good NREM period
%NREM_end  : time in seconds (from NS) of end of ~5s good NREM period
 
% output:
%newstates: a new states mat file saved with cleaner NREM 
tic
Rs = 1250;
state_mat = dir('*-states*');
load(state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5}; 
NREM = StateIntervals{3};
WAKE = StateIntervals{1};
state = NREM; %Set state as NREM

%Load CH
lfp_file = dir('*.lfp');
[~, fbasename, ~] = fileparts(lfp_file.name);
CH_Nall=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));

disp(['Cleaning NREM for ', fbasename])
% Read LFP from one channel
lfp= readmulti(lfp_file.name,CH_Nall,goodch);     
disp('Done reading LFP')
 
% Get the interval of the states
state_num = 3;
state_interval=states2interval(states,state_num)*Rs;

% Extract LFP in the NREM state in noisy band
n=3;
Wn=[200 400]; %[110 180]; %% [110 180]; %% [85 150]
[b,a]=butter(n,2*Wn/Rs,'bandpass'); 

lfp_fil= filtfilt(b,a,lfp);
[fil_restricted,fil_restricted_time]=extract_data(lfp_fil,state_interval);

%Squared filtered signal
sig = fil_restricted.^2;  

%Good NREM
ti = NREM_start*Rs;
tf = NREM_end*Rs;

tis = find(fil_restricted_time == round(ti));
tif = find(fil_restricted_time == round(tf));

%Determine std value from good NREM section
thres = std(sig(tis:tif));

%Find threshold 25*std 
noise_factor = 25;
noise= noise_factor*thres;
threshold_index=sig>noise ;

noise_time = fil_restricted_time(threshold_index); 
newstates = states;

%Change the states 
for i = 1:length(noise_time)

   res = round(noise_time(i)/1250);
   newstates(res) = 1;

end

filename= strcat(fbasename, '-newstates');
save(filename, 'newstates');


A = length(find(states==3));
B =length(find(newstates==3));
disp('NREM cleaned :)') 

disp(['Old NREM: ', num2str(A), ' new NREM: ', num2str(B)])
toc    
end

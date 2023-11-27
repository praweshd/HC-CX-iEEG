function [res_ref voltage]=Peak_finder_negative_voltage(LFP,low_lim, high_lim,Fs)    
% Usage: finds negative peaks in LFP
% PD, DKH, AH
    %% LFP=  single channel data, low_limit, high_lim=lower and higher thershold (absolute value) 

    %LFP_limit= median(abs(LFP))/0.6745;     %calcuate the noise level
    S=diff(LFP);                            % Calculate the slopes
    S1=S(2:end);
    S2=S(1:end-1);

    imn=find(S1.*S2 <= 0 & S1-S2 > 0 & S1 > 0)+1;    % find the locations of negative peaks
    imn_values=LFP(imn) ;                            % store the values of the peak
    imn_trig=  (imn_values < -1*low_lim  & imn_values> -1*high_lim); %choose the peaks matching the bandwidth
    immn=imn(find(imn_trig==1));                     % spike times
    
    res=immn; % all peak time

    %% applying delay
    
    if length(res)>2
        refreactory= 20; % avoid detection of complex waveform peaks; unit is no of samples 
%         fprintf('\nAssuming %f s refractory period\n',refreactory/Fs);
        res_diff= diff(res);
        res_refract=res_diff>refreactory; % Contains sample number of spikes that obey the refractory period
        res_tmp=res(2:end).*res_refract; % Eliminating the ones that do not obey refractory period

        %% extraction
        res_ref= res_tmp(find(res_tmp~=0)); % spikes peak time
        voltage=LFP(res_ref);
    else
        res_ref=[];
        voltage=[];
    end
    

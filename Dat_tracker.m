function [data]=Dat_tracker(F_in, T, N_points,N_channels)
% USAGE 
% Loads data from .DAT of .FIL based on given time. It's based on
% file-functions of matlab so it only opens and read and does not load, very
% convenient for large data files.
%
% INPUTS
%       f_in:          filename with extention. for example: 'data.dat'
%       T:             the time stamp (sample number) of intrest the data will be centred to time T (can be a 1D matrix)
%       N_points:      Total number of points of the waveform
%       N_channels:    number of channels in the data file 

% OUTPUT:
%       DATA:           a 3 dimnensional matrix [channel_number, data, T]

% Prawesh Dahal    

if rem(N_points,2)~=0
    shift=ceil(N_points/2);
else
    shift=(N_points/2);
end
    T= T - shift; 

if T<0
    disp('Error: not enough data');
    data=[];
else
    for i=1:length (T)
        fid=fopen(F_in,'r');  % open the file and get the ID
        pointer=ftell(fid);                               %read the current position in the file 
        pointer=fseek(fid,2*N_channels*(T(i)),'bof');    % move to the desired posistion, the file formate is : Ch1_sample1, CH2_sample1,.... ChN_sampleN
        
        if pointer==-1
            disp('Error: not a correct pointer');
            data=[];
        else
        pointer=ftell(fid);                               %check the currrent posistion in the file 

        data(:,:,i)= fread(fid,[N_channels,N_points ],'int16');   % reads the .dat file and put it into a matrix of [channel_number  data] 
        fclose(fid);
        end
        
    end
end 

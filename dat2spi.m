function [S] = dat2spi(Data_dat,resample_factor,arburg_n,Rs,freq,smooth_win)
% USAGE 
% To identify segments of EEG (.dat) during which the ratio of power in SPI band (9-20) than other high/low
% bands is high. Use this .spi file to then later detect SPI start and end time points
%
% INPUTS
%       Data_dat:      filename
%       resample_factor: new Rs              
%       arburg_n:      2
%       freq:    f band to perform wavelet (1-100)  
%       smooth_win: window in ms for gaussian smoothing

% OUTPUT:
%       S:       a 2D matrix, for each channel (row), columns are SPI ratio
%       in time

% Prawesh Dahal (revised for Brain paper 2019)  

Data_spi=   resample(Data_dat',1,resample_factor)';
[c,~]= size(Data_spi);

for i = 1: c
    x=Data_spi(i,:);
    b = arburg(x,arburg_n);
    y = Filter0In(b, x);
    [S_tmp,~,~] = awt_freqlist(y,Rs/resample_factor,freq,'Gabor');

    S_tmp_smooth = smoothdata(abs(S_tmp),1,'gaussian',smooth_win);

    spi_band = sum(S_tmp_smooth(:,9:20), 2);
    low_band = sum(S_tmp_smooth(:, 1:5), 2);
%     high_band = sum(S_tmp_smooth(:, 25:40), 2);
%     spi_ratio= (spi_band - low_band - high_band) ./ (spi_band + low_band + high_band);
    
    spi_ratio= (spi_band - low_band) ./ (spi_band + low_band);


    S(i,:)  = spi_ratio;
end
%     x=reshape(Data_spi',1,c*r);
%     b = arburg(x,arburg_n);
%     y = Filter0In(b, x);
%     [S_tmp,~,~] = awt_freqlist(y,Rs/resample_factor,freq,'Gabor');
% 
%     S_tmp_smooth = smoothdata(abs(S_tmp),1,'gaussian',smooth_win);
% 
%     spi_band = sum(S_tmp_smooth(:,10:20), 2);
%     low_band = sum(S_tmp_smooth(:, 2:8), 2);
%     high_band = sum(S_tmp_smooth(:, 25:40), 2);
%     spi_ratio= (spi_band - low_band - high_band) ./ (spi_band + low_band + high_band);
% 
% 
%     S  = reshape (spi_ratio',r,c)';
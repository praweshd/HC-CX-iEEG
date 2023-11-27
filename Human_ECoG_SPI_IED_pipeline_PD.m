% Sample pipeline for Human IED data
% Prawesh Dahal (2018)
clc
clear
close all
%% session dependent variables 
IED_name='NY394_dayfile_DAY3_AM119NREMIED1.mat';
basename='NY394_dayfile_DAY3_AM';

%% Load the grids and strips CH keys
Rs=1250;
freq_band=1:0.5:50;
load (cat(2,basename,'_CH_key.mat'));
CH_N=length(CH_key);
load NY394_map_Grid_strips
Gridmap=NY394.grid;

%% Loading IED event times  
load (IED_name);
res_IED=round(IED(:,1).*Rs);
%% Calulating mean and STD for zscoring 
load (cat(2,basename,'-states.mat'))
%%
NREM=find(states==3);
time=NREM(500)*Rs;
duration =  600*Rs;
%%
Data_zscore_all=Dat_tracker(cat(2,basename,'.lfp'),time,duration,CH_N);
Data_zscore=Data_zscore_all(1:126,:);
%% Perform wavelet in 1-50 Hz band
clear Data_zscore_all;
for CH=1:126
    [S,freq]=awt_freqlist_Frank(Data_zscore(CH,:),Rs,freq_band,'Gabor');
    pwr_spectrum.mean(CH,:)=mean(abs(S));
    pwr_spectrum.std(CH,:)=std(abs(S));
    clc
    disp(cat(2,'zscoring channel #',num2str(CH)));
end
save('pwr_spectrum','pwr_spectrum');
%% Perform SPI band filtering and hilbert
fc= [9 18];
filter_type='bandpass';
n=1;
Wn=fc;
[b,a]=butter(n,2*Wn/Rs,filter_type);
Data_fil=filtfilt(b,a,Data_zscore')';
spi_baseline_pwr.mean=  mean(   abs(hilbert(Data_fil')'));
spi_baseline_pwr.std=  std (   abs(hilbert(Data_fil')'));
save('spi_baseline_pwr','spi_baseline_pwr');

clear time 
clear duration
clear Data_fil

%% Load LFP
time=res_IED(1:end);
duration=6.*Rs;
Data= Dat_tracker(cat(2,basename,'.lfp'),time,duration,CH_N);
%% filtering and Hilbert
clear S
clear Data_zscore

fc= [9 18];
filter_type='bandpass';
n=1;
Wn=fc;
[b,a]=butter(n,2*Wn/Rs,filter_type);
for CH=1:126
   Data_fil=filtfilt(b,a,squeeze(Data(CH,:,:)));
   Data_hilbert=abs(hilbert(Data_fil));
   Data_hilbert_mean(CH,:)=mean(Data_hilbert,2);
   Data_hilbert_std(CH,:)=std(Data_hilbert,[],2)./length(res_IED);
    clc
   disp(cat(2,'filtering and hilberting channel #',num2str(CH)));
   clear Data_fill
end
%%
Spi_pwr_hilbert.mean=Data_hilbert_mean;
Spi_pwr_hilbert.std=Data_hilbert_std;
save('Spi_pwr_hilbert','Spi_pwr_hilbert')

% [channel_N data_length trial_length]=size(Data_hilbert);
% Spi_pwr_hilbert.mean=mean(abs(Data_hilbert),3);
% Spi_pwr_hilbert.std=std(abs(Data_hilbert),[],3)./sqrt(trial_length);
% save('Spi_pwr_hilbert','Spi_pwr_hilbert')
% clear Data_fil
% clear Data_hilbert
%% Time-frequency averaging 

spi_mean= pwr_spectrum.mean;
spi_std= pwr_spectrum.std;
%%
[channel_N t_total trial_N ]=size(Data);
for CH=1:126
        data_CH=squeeze(Data(CH,:,:));
        Z=ones(t_total,length(freq_band),trial_N); 
        for trial=1:trial_N;
            [S,freq]=awt_freqlist_Frank(data_CH(:,trial),Rs,freq_band,'Gabor');
             S_abs=abs(S);
             S_zscore=ones(t_total,length(freq_band)); 
             for f=1:length(freq_band)
                 for t=1:t_total
                     S_zscore(t,f)=  (  S_abs(t,f)-spi_mean(CH,f)  ) ./  spi_std(CH,f);
                 end
             end
             Z(:,:,trial)= S_zscore;
        end
            TF_avg(CH).S_zscore= mean(Z,3);
            disp(cat(2,'wavelet of channel ' , num2str((CH)), ' is done!'));
end


save ('TF_avg','TF_avg','-v7.3')


%% Plots 

power_avg=Spi_pwr_hilbert.mean;
power_std=(Spi_pwr_hilbert.std)./sqrt(648);
%% GRID 
figure_grid_spi_hilbert=figure_ctrl('Spi Hilbert @ grid',1500,1000);
for Ch=1:64
     if isnan(Gridmap(Ch)) ~= 1; 
        subaxis(8,8,Ch, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.05);
        plot(power_avg(Gridmap(Ch),:) ,'k','LineWidth',2);
        hold on 
        plot(power_std(Gridmap(Ch),:)+power_avg(Gridmap(Ch),:) ,'--r');
        plot(power_avg(Gridmap(Ch),:)-power_std(Gridmap(Ch),:) ,'--r');
        xlim ([500 7000]);
     end
end
print(figure_grid_spi_hilbert,'-djpeg ',cat(2,'figures/',figure_grid_spi_hilbert.Name));
%% strips
strip_names=fieldnames(NY394);
strip_N=length(strip_names);
%%
for strips=1:strip_N-1
    
    figure_strip_spi_hilbert=figure_ctrl(cat(2,'Spi Hilbert @',cell2mat(strip_names(strips)) ),200,1000);
    strip_channels=getfield(NY394,strip_names{strips});
    for Ch=1: length(strip_channels)
        subaxis(6,1,Ch, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.05);
        plot(power_avg(strip_channels(Ch),:) ,'k','Linewidth',3);
        hold on 
        plot(power_std(strip_channels(Ch),:)+power_avg(strip_channels(Ch),:) ,'--r');
        plot(power_avg(strip_channels(Ch),:)-power_std(strip_channels(Ch),:) ,'--r');
        xlim([500 7000])
        %axis_cleaner
    end
    print(figure_strip_spi_hilbert,'-djpeg ',cat(2,'figures/',figure_strip_spi_hilbert.Name));
    
end
close all

%%
load TF_avg.mat
%% TF of Grid 
figure_TF_IED_trg_spi_grid=figure_ctrl('TF_IED_trg_spi@ grid',2000,1000);
spi_pwr=Spi_pwr_hilbert.mean;

for Ch=1:64
    if isnan(Gridmap(Ch)) ~= 1; 
        subaxis(8,8,Ch, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
        qq=TF_avg(Gridmap(Ch)).S_zscore;
        qq_restricted=qq(1000:6500,1:50);
        imagesc(flipud(qq_restricted)' ); axis xy; caxis([0 0.2]);
        hold on 
        plot(spi_pwr(Gridmap(Ch),1000:6500),'w','LineWidth',2);
        axis_cleaner;
    end
end
 %print(figure_TF_IED_trg_spi_grid,'-djpeg ',cat(2,'figures/',figure_TF_IED_trg_spi_grid.Name));
%% TF of strips
for strips=1:strip_N-1
    
    figure_strip_spi_TF=figure_ctrl(cat(2,'TF_IED_trg_spi@',cell2mat(strip_names(strips)) ),200,1000);
    strip_channels=getfield(NY394,strip_names{strips});
    for Ch=1: length(strip_channels)
        subaxis(8,1,Ch, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.05);
        qq=TF_avg(strip_channels(Ch)).S_zscore;
        qq_restricted=qq(1000:6500,1:50);
        imagesc((qq_restricted)' ); axis xy; caxis([0 0.2]);
        hold on 
        plot(spi_pwr(strip_channels(Ch),1000:6500),'w','LineWidth',2);
        axis_cleaner;
    end
    print(figure_strip_spi_TF,'-djpeg ',cat(2,'figures/',figure_strip_spi_TF.Name));
    
end

 
%% IED LFP trigger avg
IED_avg=mean(Data,3);

%% IED avg of Grid
figure_IED_avg=figure_ctrl('IED triggered average@grid',1500,800);
for Ch=1:64
    if isnan(Gridmap(Ch)) ~= 1; 
        subaxis(8,8,Ch, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); 
        plot(IED_avg(Gridmap(Ch),:),'k','LineWidth',1);
        ylim([-40 40])
        axis_cleaner;
        axis_cleaner;
    end
end
 print(figure_IED_avg,'-djpeg ',cat(2,'figures/',figure_IED_avg.Name));
%% TF of strips
for strips=1:strip_N-1
    figure_strip_IED_avg=figure_ctrl(cat(2,'IED triggered average@',cell2mat(strip_names(strips)) ),200,1000);
    strip_channels=getfield(NY394,strip_names{strips});
    for Ch=1: length(strip_channels)
        subaxis(8,1,Ch, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.05);
        plot( IED_avg(strip_channels(Ch),:),'k','LineWidth',1); axis xy; caxis([0 0.2]);
    end
    print(figure_strip_IED_avg,'-djpeg ',cat(2,'figures/',figure_strip_IED_avg.Name));
    
end
%% CCG calculation 
load NY394_spindles_day3_PM.mat
%%
Win=10000*2;
Bin_ms=200;
res_IED=IED(:,1).*Rs;
%%
 for CH=1:126 
    spi_time= NY394_spindles(CH).spindles;
     %res_spi=round((spi_time(:,1).*Rs + spi_time(:,3).*Rs)/2);
     res_spi=spi_time(:,2).*Rs;

    if ~isempty(res_spi) 
        [H(CH,:),B ]= CrossCorr_ms( res_IED/Rs*1e3,res_spi/Rs*1e3   ,Bin_ms  ,Win/Bin_ms);
        %H_mean(CH,:)=mean(H(CH,:));
        cchEvt = (H(CH,:).*length(res_IED)*Bin_ms/1000)';
        window = 21;                %number of bins for convolution
        alpha = 0.05;

        [dumy, pred, dumy ] = cch_conv(round(cchEvt),window);
        hiBound = poissinv( 1 - alpha, pred);
        loBound = poissinv( alpha, pred);

        hiB(CH,:) = hiBound/(length(res_IED)*Bin_ms/1000);
        loB(CH,:) = loBound/(length(res_IED)*Bin_ms/1000);
    else
     H(CH,:)=zeros(21,1);
    end         
 end


 %% plot CCG of IED and spindles
 figure_ctrl('CCG_IED_SPI',1000,1000);
 for CH=1:64
    if isnan(Gridmap(CH)) ~= 1;
        subaxis(8,8,CH, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.01);
        bar(B,H(Gridmap(CH),:));
        hold on 
        plot(B,hiB(Gridmap(CH),:),'r--')
        plot(B,loB(Gridmap(CH),:),'r--')
        axis([-1e4/2 1e4/2 0 0.2])
    end

end

%%

%% CCG of strips
for strips=1:strip_N-1
    figure_strip_IED_avg=figure_ctrl(cat(2,'IED triggered average@',cell2mat(strip_names(strips)) ),100,1000);
    strip_channels=getfield(NY394,strip_names{strips});
    for Ch=1: length(strip_channels)
        %subaxis(8,1,Ch, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.05);
        bar(B,H(strip_channels(Ch),:));
        hold on 
        plot(B,hiB(strip_channels(Ch),:),'r--')
        plot(B,loB(strip_channels(Ch),:),'r--')
        axis([-1e4/2 1e4/2 0 0.3])
    end
    %print(figure_strip_IED_avg,'-djpeg ',cat(2,'figures/',figure_strip_IED_avg.Name));
    
end

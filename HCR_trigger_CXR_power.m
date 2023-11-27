%% USAGE 
% HC-ripple peak triggered power 
% Identify at the time of ripples, the power in other f brands in other brain regions
% Prawesh Dahal (May 2022) 
close all
clear all
clc 

close all; clear all; clc;
ani = 'OR26'; PPC_CH = 120;
load([ani,'_filesorder.mat']);
totalses = max([files.sesorder]);
dirs = pwd;
mkdir TriggerPower


for ses = 9
    
    ind = find([files.sesorder] == ses);
    
    sesfolder = strcat(files(ind).folder, '/', files(ind).name);
    cd(sesfolder)
    disp(['Matrix now for ', files(ind).name])  
    
    detec_file = dir('*_RIP_RES.mat'); load(detec_file.name);
    lfp_file = dir('*.lfp'); lfp_filename=lfp_file.name;
    
    detec_file = dir('*_HC_NREM_ripples.mat');
    if isempty(detec_file)
        detec_file = dir('*_HC_NREM.mat');
    end
    load(detec_file.name);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HC = ripples_Hc(:,2);
    
    CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));    
    
    Rs = 1250;
    mylfp = ripples;
    CH = PPC_CH ; %regCH(i).chs(j);
    
%     res = (mylfp(CH).res(:,2));
    res = HC;
    trgtrace = 0;    spectrace = 0;

    tottrial = length(res);%500;
    
    for j = 1:tottrial
        
        timep = res(j);
        time = round(timep*Rs);
        duration=0.2*Rs;
        dataLFPall=Dat_tracker(lfp_file.name,time,duration,CH_N);         
        dataLFP = dataLFPall(PPC_CH,:);
        
%         cxripdata(j,:) = dataLFP;
        
        arburg_n = 2;
        b = arburg(dataLFP,arburg_n);
        y = Filter0In(b, dataLFP);   
        [SS,freq,~] = awt_freqlist(y,Rs,0:250,'Gabor');
        dd = abs(SS);
        trgtrace = trgtrace + y; 
        spectrace = spectrace + dd;
    end   


spectrace = spectrace./tottrial;
xtime = -0.1:(1/1250):0.1;
xtime = xtime(1:end-1);
yy = trgtrace./tottrial;


close all
figure_ctrl('FIL',700,900);

subplot(211)  
plot(xtime,yy)
ylim([-40 40])
xlim([-0.05 0.05]) 

center = round(length(dataLFP)/2); 
around = round(0.01*Rs); 
regionpwr = spectrace(center-around:center+around, 111:181); 

subplot(212)
imagesc(xtime,0:250,spectrace'); axis xy; colormap jet;
xlim([-0.05 0.05])
caxis([1.6 2.4])
ylim([50 250])
colorbar
title([num2str(ses), '- ', files(ind).name,' ', files(ind).tags, ' = ' num2str(mean(mean(regionpwr)))],'Interpreter','None')

% print([files(ind).name,'_hccxtrg'],'-dpng') 
% movefile([strcat(files(ind).name,'_hccxtrg', '.png')], strcat(dirs,'/','TriggerPower'))    
% 
% print([files(ind).name,'_hccxtrg.eps'], '-depsc', '-tiff','-r300', '-painters')
% movefile([files(ind).name,'_hccxtrg.eps'], strcat(dirs,'/','TriggerPower'))

cd ..





end
%% Classify local vs global SPI based on power
% Prawesh Dahal
 
% cd '/media/prawesh/JG2/Jose/RatD21/20200803_B1'; map = 'v4';
% timep = 7269.713 ; %7269.713 %7197.833; % 7308.777; %;7269.713; %7197.833; %7130.618;
% %global 8474.707; 15460.400
% clusfile = dir('*_spiclus.mat'); load(clusfile.name);


%%%%%%%%%%%%%%--------%%%%%%%%%%%%%%%%%%%%%%%%
%OR17
close all; clear all; clc;
cd '/media/prawesh/73c98bba-4822-49b8-8979-c608a3adee84/RAT/ORat_17/20200220_0um_posttrain'
% local timep's - 6951.102; 6877.645*R*; 6886.144*S*; 6888.777*S*; 7174.215;
% global timep's - 7389.860; 2844.981; 4646.684* ; 5598.023;;
% ? 5760.622; 6039.500
timep = 11993.7708 ; %11993.7708; %2865.03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OR20
% clusfile = dir('*_spi_clus_nw.mat'); load(clusfile.name);
% [globidx,locidx] = glob_vs_local(spiclus,1:length(spiclus));
%%%%%%%%%%%%%%--------%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clc
% load('20200220_0um_posttrain_SPI_clus22.mat')
% clearvars -except powSheet 
% sheetX = sortrows(powSheet, 6);  
% powsheetidx = 135;
% timep = sheetX(powsheetidx,2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timep = spiclus(28).res;  
typ = 'G'; 

map = 'v4';
Rs=1250; 
duration = 1*Rs;
lfp_file=dir('*.lfp');
basename=lfp_file.name(1:end-4);
CH_N=xml2CH_N(cat(2,basename,'.xml'));

[NGmap,NGdim] = loadNGmap(map); NGlayout = reshape(NGmap,NGdim(1),NGdim(2))';

time = round(timep*Rs);
            
dataLFP=Dat_tracker(lfp_file.name,time,duration,CH_N);

NGcol = 5; 
CH_sel_all = NGlayout(:,NGcol);
CH_sel_NS = CH_sel_all; 

%For OR17
% CH_sel_NS = [45 48 8 98 24 80 119 90 50 47 35 110 21 7 150 153] +1; %this was for ant1

% CH_sel_NS = [48 38 68 82 26 78 18 96 47 127 99 6 93 66] +1; %this was for ant1

% CH_sel_NS = [45 48 95 64 13 9 90 40 54 50 35 1 28 21] +1; %this was for rsc sheet 71

% CH_sel_NS = [45 48 97 98 20 79 40 119 54 50 2 122 66 87] +1; %this was for G1


bad_ch = [];

figure_ctrl('FilHil',2000,700);
subplot(2,5,[1 6])
%Plot the raw traces COL WISE
data=dataLFP(CH_sel_all,:); 
for ii = 1:length(CH_sel_all)
    if ~ismember(CH_sel_all(ii),bad_ch)
        plot (data(ii,:) -ii*2000,'k'); hold on
        plot (duration/2, data(ii, length(data)/2) -ii*2000, 'r.')
        axis_cleaner; axis tight
        hold on
    end
end
title([typ, ': ', num2str(timep)])

figure(1)
subplot(2,5,[3 8])
%Plot the raw traces NS WISE
data=dataLFP(CH_sel_NS,:);  
for ii = 1:length(CH_sel_NS)
    if ~ismember(CH_sel_NS(ii),bad_ch)
        plot (data(ii,:) -ii*2000,'k'); hold on
        plot (duration/2, data(ii, length(data)/2) -ii*2000, 'r.')
        axis_cleaner; axis tight
        hold on
    end
end
title([typ, ': ', num2str(timep)]) 

n=3; Wn=[10 20];
[b,a]=butter(n,2*Wn/Rs,'bandpass');
datafil = filtfilt(b,a,dataLFP')';
datahil = hilbert(datafil')';
  
subplot(2,5,[2 7])
k=0;
for ii = 1:length(CH_sel_all)
    if ~ismember(CH_sel_all(ii),bad_ch)
        k=k+1;
        plot( abs(datahil(CH_sel_all(ii),:))-k* 2000 ,'r'); 
        hold on
        plot(  datafil(CH_sel_all(ii),:) -k* 2000 ,'b'); 
        axis_cleaner; axis tight
        hold on
    end
end

%For NS
subplot(2,5,[4 9])
k=0;
for ii = 1:length(CH_sel_NS)
    if ~ismember(CH_sel_NS(ii),bad_ch)
        k=k+1;
        plot( abs(datahil(CH_sel_NS(ii),:))-k* 2000 ,'r'); 
        hold on
        plot(  datafil(CH_sel_NS(ii),:) -k* 2000 ,'b'); 
        axis_cleaner; axis tight
        hold on
    end
end
% title(['labels ', num2str(sheetX(powsheetidx,3))])
 
SPIpower = max(abs(datahil(1:128,:)),[],2); 
SPIpower2 = SPIpower(NGmap);

occ_map = reshape(SPIpower2,NGdim(1),NGdim(2))';
subplot(2,5,[5])
h= imagesc( occ_map ); hold on
set(h,'alphadata',~isnan(occ_map))
colormap jet; colorbar;
[spext_val] =  weightedSPIext(NGdim, occ_map); 
title(['column ', num2str(NGcol)])
clim([250 600])

std_file = dir('*_LG_stdminmax.mat');
load(std_file.name)
subplot(2,5,10)
plot(sort(SPIpower2,'descend'))
hold on; plot(1:length(SPIpower2), CHstd*3*ones(length(SPIpower2),1))
xlabel('Channels')
ylabel('Power')
%%

svdir = '/home/prawesh/Documents/Insync/Prawesh_PPU/Prawesh_GlobalSPI/EPS';
print(['RawTrace_March_global.eps'], '-depsc', '-tiff','-r300', '-painters')
movefile(['RawTrace_March_global.eps'], svdir)

%% Another sample plot of raw trace + wavelet + filtered + std

close all; clear all; clc;
cd '/media/prawesh/73c98bba-4822-49b8-8979-c608a3adee84/RAT/ORat_17/20200220_0um_posttrain'
 
timep = 4447.132+0.1 ; %Ask dion2828.017 -0.1;
CHsel = 111; %68; %29; %55; 
map = 'v4';
Rs=1250; 
duration = 2.5*Rs;
lfp_file=dir('*.lfp');
basename=lfp_file.name(1:end-4);
CH_N=xml2CH_N(cat(2,basename,'.xml'));

 
time = round(timep*Rs);
            
dataLFP=Dat_tracker(lfp_file.name,time,duration,CH_N);

subplot(311)
data=dataLFP(CHsel,:);
plot (data(1,:),'k'); hold on
xlim([0 duration])
ylim([-2500 2500])
title([num2str(timep), '-CH', num2str(CHsel)])

arburg_n =1; freq = 1:80;
b = arburg(data,arburg_n);
y = Filter0In(b, data);
[SS,freq,~] = awt_freqlist(y,Rs,freq,'Gabor');
dd = abs(SS);


subplot(312)
imagesc(dd' ); axis xy; colormap jet;
ylim([0 50]);  xlim([0 length(data)])
xlim([0 duration])
clim([0 20])
 
 n=3; Wn=[10 20];
[b,a]=butter(n,2*Wn/Rs,'bandpass');
datafil = filtfilt(b,a,dataLFP')';
datahil = hilbert(datafil')';

subplot(313)
plot( abs(datahil(CHsel,:))  ,'r');
hold on
plot(  datafil(CHsel,:)  ,'b'); hold on
load('20200220_0um_posttrain_LG_stdminmax.mat')
plot(1:length(data), ones(length(data),1).*CHstd, '--k')
xlim([0 duration])


plot(1:length(data), ones(length(data),1).*3*CHstd, '-g')
xticks([(0:625:duration)])
xlim([0 duration]); ylim([-1000 1000])
%%
%Normalize 
% mu = 583.6156; stD = 135.9110; 
% SPIpowerN = (SPIpower2 - mu)./stD;

% pmin = 325.6486 ; pmax = 1.3390e+03;
% SPIpowerN = (SPIpower2 - pmin)./(pmax - pmin);

% pmin = 144.8328; pmax = 1.8872e+03;
% SPIpowerN = (SPIpower2 - pmin)./(pmax - pmin);
% 
% occ_mapN = reshape(SPIpowerN,NGdim(1),NGdim(2))';
% subplot(2,4,[7])
% h= imagesc( occ_mapN ); hold on
% set(h,'alphadata',~isnan(occ_mapN))
% colormap jet; colorbar;
% [spext_valN] =  weightedSPIext(NGdim, occ_mapN); 
% title(num2str(spext_valN))
% 
% subplot(2,4,4)
% histogram(SPIpower2,6) 
% 
% subplot(2,4,8)
% histogram(SPIpowerN,6) 

% Electrode spatial extent dist-weighted sum

%% Try this weighted-dist for all SPI 

close all
clear all
clc

clusfile = dir('*_spi_clus_nw.mat'); load(clusfile.name); 

map = 'v4'; Rs=1250; duration = 1*Rs; 
lfp_file=dir('*.lfp');
basename=lfp_file.name(1:end-4);
CH_N=xml2CH_N(cat(2,basename,'.xml'));

[NGmap,NGdim] = loadNGmap(map); NGlayout = reshape(NGmap,NGdim(1),NGdim(2))'; 
k = 0; 

spiclus = spiclus(2:end);
[minP, maxP] = minmaxpow(spiclus, NGmap);  

%Form labels
labels = zeros(length(spiclus),1);
thresh = 0.55;

for niter = 1:length(spiclus) 
     
    timep = spiclus(niter).res ; 
    time = round(timep*Rs);
    
    dataLFP=Dat_tracker(lfp_file.name,time,duration,CH_N);
    
    n=3; Wn=[10 20];
    [b,a]=butter(n,2*Wn/Rs,'bandpass');
    datafil = filtfilt(b,a,dataLFP')';
    datahil = hilbert(datafil')';
    
    SPIpower = max(abs(datahil(1:128,:)),[],2);
    SPIpower2 = SPIpower(NGmap);

    %Normalize
    %mu = 583.6156; stD = 135.9110;
    %SPIpowerN = (SPIpower2 - mu)./stD;
%     pmin = 325.6486 ; pmax = 1.3390e+03;
%     pmin = 144.8328; pmax = 1.8872e+03;
    SPIpowerN = (SPIpower2 - minP)./(maxP - minP); 
    
    k = k + 1;
    pcaX(k,:) = SPIpowerN;
    occ_map = reshape(SPIpowerN,NGdim(1),NGdim(2))';
    
    [spext_val] =  weightedSPIext(NGdim, occ_map);
     
    val(niter,1) = spext_val;
    val(niter,2) = min(SPIpower2);
    val(niter,3) = max(SPIpower2);
     
    %Label Method 1 - formula
%     if spext_val > thresh
%         labels(niter) = 1; %it becomes global
%     end
%     
    %Label Method 2 - Detection Power
    stdLim = 150.7023 *3;
    crit = find(SPIpower2 < 465);
    if isempty(crit) 
        labels(niter) = 1; %no low power value call it global
%     else
%         if length(crit) < 6
%             labels(niter) = 1; 
%         else 
%             labels(niter) = 0;
%         end
    end
    
    powSheet(niter,1) = niter ; 
    powSheet(niter,2) = labels(niter);
    powSheet(niter,3) = min(SPIpower2);
    powSheet(niter,4) = max(SPIpower2); 
    powSheet(niter,5) = spext_val;  
     
end

 
%% Do PCA  

X = pcaX;
[coeff,score,latent,tsquared,explained,mu] = pca(X);

close all
figure_ctrl('PC', 600, 500)
hold on
subplot(131)
scatter(score(:,1),score(:,2),10,'k','filled');

subplot(132)
hold on
scatter3(score(labels==0,1), score(labels==0,2), score(labels==0,3),10,'b','filled')
scatter3(score(labels==1,1), score(labels==1,2), score(labels==1,3),10,'r','filled')
legend('Local','Global')
xlabel('PCA1'); ylabel('PCA2'); box on;

detec_file = dir('*_HC_NREM_ripples.mat'); load(detec_file.name); hcres = ripples_Hc(:,2); 
spires = [spiclus.res];

[~,~, ~, ~, ~,matchedidx] = LFPmatching(hcres,spires,300);

%%
subplot(133)
hold on
scatter(score(:,1),score(:,2),10,'k','filled'); 
scatter(score(matchedidx,1),score(matchedidx,2),10,'g','filled'); 
legend('HC-unMatched','HC-Matched')
xlabel('PCA1'); ylabel('PCA2'); box on;




%% Try density plot
close all

H = densityplot(score(:,1),score(:,2));

figure(2)
plot(score(:,1),score(:,2),'k*','MarkerSize',5);

 
%% try Kmeans

KmeansData = [score(:,1) score(:,2)];
close all
figure;
plot(KmeansData(:,1),KmeansData(:,2),'k*','MarkerSize',5);

opts = statset('Display','final');
[idx,C] = kmeans(KmeansData,2,'Distance','cityblock',...
    'Replicates',5,'Options',opts);

figure;
plot(KmeansData(idx==1,1),KmeansData(idx==1,2),'r.','MarkerSize',12)
hold on
plot(KmeansData(idx==2,1),KmeansData(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off


%%
close all
figure_ctrl('Coeff', 1800, 900)
 for i = 1:10
    subplot(2,5,i)
    pcamap =   reshape(coeff(:,i),NGdim(1),NGdim(2))';
    h= imagesc( pcamap ); hold on
    set(h,'alphadata',~isnan(pcamap))
    colormap jet; colorbar;
end

%%
close all
figure_ctrl('Obs',800,600)
for i = 1:3
    subplot(2,2,i)
    hold on
    %     histogram(score(:,i))
    histogram(score(labels==0,i),'Normalization','Probability','BinWidth',0.2 )
    histogram(score(labels==1,i),'Normalization','Probability','BinWidth',0.2 )
    %     legend('Local', 'Global')
    
    title(['Var: ', num2str(explained(i))])
    

end

%%
close all
figure_ctrl('PC', 600, 500)
hold on
% scatter3(score(labels==0,1), score(labels==0,2), score(labels==0,3),10,val(labels==0,1),'filled')
% scatter3(score(labels==1,1), score(labels==1,2), score(labels==1,3),10,val(labels==1,1),'filled')
% colormap(jet)
% % set(gca,'ColorScale','log')
% colorbar

scatter3(score(labels==0,1), score(labels==0,2), score(labels==0,3),10,'b','filled')
scatter3(score(labels==1,1), score(labels==1,2), score(labels==1,3),10,'r','filled')
legend('Local','Global')
% title(num2str(thresh))
xlabel('PCA1'); ylabel('PCA2'); box on;
%% tsne
Y = tsne(X,'Algorithm','exact','NumPCAComponents',5,'Distance','chebychev');
close all
figure()
gscatter(Y(:,1),Y(:,2),labels);

%%
close all
figure_ctrl('FilHil',1800,800);
subplot(121)
histogram(val(:,1),100)  
% subplot(132)
% histogram(val(:,2),60) 
% ylabel('Counts'); xlabel('Mean SPI power per event');
% title(['Mean SPI Power, \mu = ',num2str(mean(val(:,2))), ' \sigma = ', num2str(std(val(:,2)))])
% 
% subplot(133)
% histogram(val(:,3),100) 
%% Predict trained model on test 
close all; clear all; clc;
cd '/media/prawesh/73c98bba-4822-49b8-8979-c608a3adee84/RAT/ORat_17/20200220_0um_posttrain'
load('trainedModel_1.mat');

cd '/media/prawesh/73c98bba-4822-49b8-8979-c608a3adee84/RAT/ORat_17/20200218_0um_posttrain'
load('LG_X_Test.mat');
load('LG_Y_Test.mat');

yfit = trainedModel020.predictFcn(X);
C = confusionmat(labels,yfit);
confusionchart(C)

%% Is the test results correct
clearvars -except powSheet labels powsheet2
close all; clc;
clusfile = dir('*_spi_clus_nw.mat'); load(clusfile.name);

map = 'v4'; Rs=1250; duration = 1*Rs;
lfp_file=dir('*.lfp');
basename=lfp_file.name(1:end-4);
CH_N=xml2CH_N(cat(2,basename,'.xml'));

[NGmap,NGdim] = loadNGmap(map); NGlayout = reshape(NGmap,NGdim(1),NGdim(2))';
k = 0;

spiclus = spiclus(2:end);
[minP, maxP] = minmaxpow(spiclus, NGmap);

for niter = 1554
    
    timep = spiclus(niter).res ;
    time = round(timep*Rs);
    
    dataLFP=Dat_tracker(lfp_file.name,time,duration,CH_N);
    
    NGcol = 8;
    CH_sel_all = NGlayout(:,NGcol);
    data=dataLFP(CH_sel_all,:);
    
    bad_ch = [];
    
    figure_ctrl('FilHil',1800,800);
    subplot(2,4,[1 5])
    %Plot the raw traces
    
    for ii = 1:length(CH_sel_all)
        if ~ismember(CH_sel_all(ii),bad_ch)
            plot (data(ii,:) -ii*2000,'k'); hold on
            plot (duration/2, data(ii, length(data)/2) -ii*2000, 'r.');
            axis_cleaner; axis tight
            hold on
        end
    end
    title([': ', num2str(timep)])
    
    n=3; Wn=[10 20];
    [b,a]=butter(n,2*Wn/Rs,'bandpass');
    datafil = filtfilt(b,a,dataLFP')';
    datahil = hilbert(datafil')';
    
    load('20200220_0um_posttrain_LG_stdminmax.mat')

    
    subplot(2,4,[2 6])
    k=0;
    for ii = 1:length(CH_sel_all)
        if ~ismember(CH_sel_all(ii),bad_ch)
            k=k+1;
            plot( abs(datahil(CH_sel_all(ii),:))-k* 2000 ,'r');
            hold on
            plot(  datafil(CH_sel_all(ii),:) -k* 2000 ,'b'); hold on
            plot((CHstd*3).*(ones(1250,1)) -k*2000, 'g--')
            axis_cleaner; axis tight
            hold on
        end
    end
    
    SPIpower = max(abs(datahil(1:128,:)),[],2);
    SPIpower2 = SPIpower(NGmap);
    SPIpowerN = (SPIpower2 - minP)./(maxP - minP);
    
    occ_map = reshape(SPIpower2,NGdim(1),NGdim(2))';
    
    subplot(2,4,[3])
    h= imagesc( occ_map ); hold on
    set(h,'alphadata',~isnan(occ_map))
    colormap jet; colorbar;
    [spext_val] =  weightedSPIext(NGdim, occ_map);
    title(['column ', num2str(NGcol)])
    
    
    occ_mapN = reshape(SPIpowerN,NGdim(1),NGdim(2))';
    subplot(2,4,[7])
    h= imagesc( occ_mapN ); hold on
    set(h,'alphadata',~isnan(occ_mapN))
    colormap jet; colorbar;
    [spext_valN] =  weightedSPIext(NGdim, occ_mapN);
    title(num2str(spext_valN))
    
    subplot(2,4,4)
    histogram(SPIpower2,6)
    title(['Prediction: ' num2str(labels(niter))])
    
    subplot(2,4,8)
    histogram(SPIpowerN,6) 
    title(['spiclus #: ' num2str(niter)])
    min(SPIpower2)
    max(SPIpower2)
end
%%


%%
close all
clusfile = dir('*_spi_clus_nw.mat'); load(clusfile.name); 
detec_file = dir('*_HC_NREM_ripples.mat'); load(detec_file.name); hcres = ripples_Hc(:,2); 
spires = [spiclus.res];
    
%Extract global/local spindles
window=20000; Rs=1250; bin=300; tol = 600;
[globidx,locidx] = glob_vs_local(spiclus,[1:length(spiclus)]);
globspi = spires(labels==1);
locspi = spires(labels==0);

% globspi = spires(globidx);
% locspi = spires(locidx);

res1 = hcres; res2 = locspi;
[H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs);
H = smooth(H,3);
[mod1] = CCGmodul(H,B, hiB, loB,tol,'midd');
figure_ctrl('A',1800,500)
subplot(131)
bar(B,H);
% title(['HCR-SPI ',num2str(length(hcres)),'-',num2str(length(PPCSPI))])
title(num2str(mod1))

hold on; plot(B,hiB,'--r','LineWidth',0.5);
hold on; plot(B,loB,'--r','LineWidth',0.5);
hold on; plot(B,(hiB+loB)./2,'--g','LineWidth',0.5);  ylim([0 0.6])
clear res1 res2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res1 = hcres; res2 = globspi;
[H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs);
H = smooth(H,3);
[mod2] = CCGmodul(H,B, hiB, loB,tol,'midd');

subplot(132)
bar(B,H); 
% title(['CXR-SPI ',num2str(length(PPCCXR)),'-',num2str(length(PPCSPI))])
title(num2str(mod2))

hold on; plot(B,hiB,'--r','LineWidth',0.5);
hold on; plot(B,loB,'--r','LineWidth',0.5);
hold on; plot(B,(hiB+loB)./2,'--g','LineWidth',0.5);  ylim([0 0.6])
clear res1 res2

subplot(133)
plot([1 2], [mod1 mod2], '-r'); hold on 

[matchedrefL,~, ~, ~, ~,~] = LFPmatching(hcres,locspi,300);
mod3 = length(matchedrefL)/length(locspi);

[matchedrefG,~, ~, ~, ~,~] = LFPmatching(hcres,globspi,300);
mod4 = length(matchedrefG)/length(globspi);

ylabel('Coupling Strength(a.u.)')

plot([1 2], [mod3 mod4], '-b'); hold on
 
xlim([0 3]); xticks([1 2]); xticklabels({'L', 'G'})
legend('CCG M', 'Joint Count');
scatter([1 2], [mod1 mod2], 40,'filled','r')
scatter([1 2], [mod3 mod4], 40,'filled','b')
 ylim([0 0.4])

%% Get mean std power from detection

function CHstd = getCHstd(spindle_channel)
 

lfp_file = dir('*.lfp');  [~, fbasename, ~] = fileparts(lfp_file.name);
CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));

filename = dir('*.lfp');
[pathstr, fbasename, fileSuffix] = fileparts(filename.name);

% spindle_channel = 87; %SomatoCH 87

state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
NREM = or(StateIntervals{2}, StateIntervals{3}); 
state = NREM;

Rs=1250;
disp('Loading lfp calculating std')
lfp = LoadLfp(fbasename,CH_N,spindle_channel);  
sleep=[Range(Restrict(lfp, state), 's') Data(Restrict(lfp, state))];

signal=sleep(:,2);
n=3;
Wn=[10 20]; %[110 180]; %% [110 180]; %% [85 150]
[b,a]=butter(n,2*Wn/Rs,'bandpass'); 

fil_sleep_V=filtfilt (b,a,signal);
datahil = hilbert(fil_sleep_V);

%%
mean(abs(datahil))
std(abs(datahil))

end 
%%

%Get the min max power
function [minP, maxP] = minmaxpow(spiclus, NGmap)

map = 'v4'; Rs=1250; duration = 1*Rs;
lfp_file=dir('*.lfp');
basename=lfp_file.name(1:end-4);
CH_N=xml2CH_N(cat(2,basename,'.xml'));

for niter = 1:length(spiclus)
    
    timep = spiclus(niter).res ;
    time = round(timep*Rs);
    
    dataLFP=Dat_tracker(lfp_file.name,time,duration,CH_N);
    
    n=3; Wn=[10 20];
    [b,a]=butter(n,2*Wn/Rs,'bandpass');
    datafil = filtfilt(b,a,dataLFP')';
    datahil = hilbert(datafil')';
    
    SPIpower = max(abs(datahil(1:128,:)),[],2);
    SPIpower2 = SPIpower(NGmap);
    
    val(niter,1) = min(SPIpower2);
    val(niter,2) = max(SPIpower2);
end

minP = min(val(:,1));
maxP = max(val(:,2));

end


function [spext_val] =  weightedSPIext(NGdim, occ_map)

Gsum = 0;
cnt = 0;
for ki = 1:NGdim(2)
    for kj = 1:NGdim(1)
        
        mainpos = [ki kj];
        dist = calcdist(NGdim, mainpos);
        mainpow = occ_map(ki,kj); 
        
        for kki = 1:NGdim(2)
            for kkj = 1:NGdim(1) 
                refpow = occ_map(kki,kkj);
                Gsum = Gsum + (mainpow * refpow * dist(kki, kkj));
                cnt = cnt + 1; 
            end
        end
        
    end
end

spext_val = Gsum/cnt;
end




function [dist ] = calcdist(NGdim, mainpos)
for i = 1:NGdim(2)
    for j = 1:NGdim(1)
        dist(i,j) = norm(mainpos - [i j]);
    end
end 
end

clear

% regional spindle-spindle coupling
SPI_files = dir('*SPIonly_cat.mat');
IED_files = dir('*_IED_cat.mat');
bin_cnt = 1;
zero_bin = round(301/2);
for pp = 1:length(SPI_files)
    load(SPI_files(pp).name);
    load(IED_files(pp).name);
    
    for jj = 1:length(Spindles_only)
        if ~isempty(IED(jj).XYZ)
         Spindles_only(jj).location = IED(jj).location;
         Spindles_only(jj).XYZ = IED(jj).XYZ;
        end
        if ~isempty(Spindles_only(jj).res) && ~isempty(Spindles_only(jj).SPI_SPI_CCG_H)
            temp_H = Spindles_only(jj).SPI_SPI_CCG_H;
            [r, ~] = size(temp_H);
            for i = 1:r
                 smooth_H(i, :) = smooth(temp_H(i,:));
            end
            temp_hiB = Spindles_only(jj).SPI_SPI_CCG_hiB.*1.1;
            SPI_CCG=sum(smooth_H(:,zero_bin-bin_cnt:zero_bin+bin_cnt),2);
            SPI_CCG_hi= sum(temp_hiB(:,zero_bin-bin_cnt:zero_bin+bin_cnt),2);
            delta_SPI_CCG= SPI_CCG -SPI_CCG_hi;
            sig_corr = find(delta_SPI_CCG > 0);
                if ~isempty(IED(jj).XYZ)
                    distance_corr = [];
                    location_corr = {};
                    amount_corr = [];
                    for ii = 1:length(sig_corr)
                        if ~isempty(IED(sig_corr(ii)).XYZ)
                        location_corr{ii} = IED(sig_corr(ii)).location;
                        temp_coord1 = IED(jj).XYZ;
                        temp_coord2 = IED(sig_corr(ii)).XYZ;
                        distance_corr(ii) = sqrt((temp_coord1(1)-temp_coord2(1))^2 + (temp_coord1(2)-temp_coord2(2))^2 + (temp_coord1(3)-temp_coord2(3))^2);
                        amount_corr(ii) = delta_SPI_CCG(sig_corr(ii));
                        clear temp_coord1 temp_coord2
                        else
                            location_corr{ii} = 'unknown';
                            distance_corr(ii) = nan;
                        end
                    end
                    distance_corr(isnan(distance_corr)) = [];
                    distance_corr(distance_corr == 0) = [];
                        Spindles_only(jj).corr_loc = location_corr;
                        Spindles_only(jj).corr_dist = distance_corr;
                        Spindles_only(jj).corr_amt = amount_corr;
                        Spindles_only(jj).delta_SPI_SPI = delta_SPI_CCG;
                        clear temp_H smooth_H temp_hiB SPI_CCG SPI_CCG_hi delta_SPI_CCG distance_corr location_corr sig_corr amount_corr
                end
        clear temp_H smooth_H temp_hiB SPI_CCG SPI_CCG_hi delta_SPI_CCG distance_corr location_corr sig_corr amount_corr        
        end
     
    end
    save(SPI_files(pp).name, 'Spindles_only');
end

%% Establish grand physio and patho zones
clear
IED_files = dir('*_IED_new.mat');

for kk = 1:length(IED_files)
   load(IED_files(kk).name)
   basename = IED_files(kk).name(1:end-13);
   m = 1;
   n = 1;
   for nn = 1:length(IED)
      if IED(nn).spicoupling == 0
           IEDSPI_physio(n) = IED(nn).channel;
           n = n+1;
      elseif IED(nn).spicoupling > 0 
          IEDSPI_patho(m) = IED(nn).channel;
          m = m+1;
      end
   end
    PhysPath.physl = IEDSPI_physio;
    PhysPath.patho = IEDSPI_patho;
    save(strcat(basename, 'PhysPath'), 'PhysPath')
    clear PhysPath IEDSPI_physio IEDSPI_patho
end

%% calculate spindle extent using grand physio and patho zones
clear
SPI_files = dir('*SPIonly.mat');
IED_files = dir('*PhysPath.mat');

for pp = 1:length(SPI_files)
    load(SPI_files(pp).name);
    load(IED_files(pp).name);
    n = 1;
    m = 1;
    for jj = 1:length(Spindles_only)
        if ismember(Spindles_only(jj).channel, PhysPath.physl) == 1 && ~isempty(Spindles_only(jj).corr_dist)
            phys_num(n) = length(Spindles_only(jj).corr_dist);
            phys_dist(n) = mean(Spindles_only(jj).corr_dist);
            phys_sigcorr(n) = mean(Spindles_only(jj).corr_amt);
            phys_loc(n).all = unique(Spindles_only(jj).corr_loc);
            n = n+1;
        elseif ismember(Spindles_only(jj).channel, PhysPath.patho) == 1 && ~isempty(Spindles_only(jj).corr_dist)
            path_num(m) = length(Spindles_only(jj).corr_dist);
            path_dist(m) = mean(Spindles_only(jj).corr_dist);
            path_sigcorr(m) = mean(Spindles_only(jj).corr_amt);
            path_loc(m).all = unique(Spindles_only(jj).corr_loc);
            m = m+1;
        end
    end
    PhysPath.phnum = phys_num;
    PhysPath.phdis = phys_dist;
    PhysPath.phcorr = phys_sigcorr;
    PhysPath.phlo = phys_loc;
    PhysPath.panum = path_num;
    PhysPath.padis = path_dist;
    PhysPath.pacorr = path_sigcorr;
    PhysPath.palo = path_loc;
    save(IED_files(pp).name, 'PhysPath');
    clear phys_num phys_dist phys_sigcorr phys_loc path_num path_dist path_sigcorr path_loc
end

%% PhysPath stats
clear
IED_files = dir('*PhysPath.mat');
for tt = 1:length(IED_files)
    load(IED_files(tt).name)
    basename = IED_files(tt).name(1:end-12);
    
    Stats.phys.num.mean = mean(PhysPath.phnum);
    Stats.phys.num.median = median(PhysPath.phnum);
    Stats.phys.num.std = std(PhysPath.phnum);
    Stats.phys.num.ste = ste(PhysPath.phnum);

    Stats.phys.dist.mean = mean(PhysPath.phdis);
    Stats.phys.dist.median = median(PhysPath.phdis);
    Stats.phys.dist.std = std(PhysPath.phdis);
    Stats.phys.dist.ste = ste(PhysPath.phdis);

    Stats.phys.corr.mean = mean(PhysPath.phcorr);
    Stats.phys.corr.median = median(PhysPath.phcorr);
    Stats.phys.corr.std = std(PhysPath.phcorr);
    Stats.phys.corr.ste = ste(PhysPath.phcorr);

    temp = PhysPath.phlo;
    list_all = {};
    for xx = 1:length(temp)
        num(xx) = length(temp(xx).all);
        list_all = horzcat(list_all, temp(xx).all);
    end
    Stats.phys.loc.mean = mean(num);
    Stats.phys.loc.median = median(num);
    Stats.phys.loc.std = std(num);
    Stats.phys.loc.ste = ste(num);
    Stats.phys.loc.uni = unique(list_all);
    PhysPath.phlonum = num;
    
    clear temp list_all num

    Stats.path.num.mean = mean(PhysPath.panum);
    Stats.path.num.median = median(PhysPath.panum);
    Stats.path.num.std = std(PhysPath.panum);
    Stats.path.num.ste = ste(PhysPath.panum);

    Stats.path.dist.mean = mean(PhysPath.padis);
    Stats.path.dist.median = median(PhysPath.padis);
    Stats.path.dist.std = std(PhysPath.padis);
    Stats.path.dist.ste = ste(PhysPath.padis);

    Stats.path.corr.mean = mean(PhysPath.pacorr);
    Stats.path.corr.median = median(PhysPath.pacorr);
    Stats.path.corr.std = std(PhysPath.pacorr);
    Stats.path.corr.ste = ste(PhysPath.pacorr);

    temp = PhysPath.palo;
    list_all = {};
    for xx = 1:length(temp)
        num(xx) = length(temp(xx).all);
        list_all = horzcat(list_all, temp(xx).all);
    end
    Stats.path.loc.mean = mean(num);
    Stats.path.loc.median = median(num);
    Stats.path.loc.std = std(num);
    Stats.path.loc.ste = ste(num);
    Stats.path.loc.uni = unique(list_all);
    PhysPath.palonum = num;
    
    clear temp list_all num

stat_name = strcat(basename, '_PhysPath_stats');
save(stat_name, 'Stats');

 % Statistical testing - nonparametric t-test
 % 1 = uncoupled
 % 2 = coupled

 [p,h,stats] = ranksum(PhysPath.phnum, PhysPath.panum);
 Stats.ttest.num.p = p;
 Stats.ttest.num.stats = stats;

 [p,h,stats] = ranksum(PhysPath.phdis, PhysPath.padis);
 Stats.ttest.dis.p = p;
 Stats.ttest.dis.stats = stats;

 [p,h,stats] = ranksum(PhysPath.phcorr, PhysPath.pacorr);
 Stats.ttest.corr.p = p;
 Stats.ttest.corr.stats = stats;

 [p,h,stats] = ranksum(PhysPath.phlonum, PhysPath.palonum);
 Stats.ttest.loc.p = p;
 Stats.ttest.loc.stats = stats;

save(stat_name, 'Stats');
clear Stats
end

%% SPI_SPI coupling segregrated by Physo Patho zones

% PhysoPatho_files = dir('*PhysPath.mat');
% load(PhysoPatho_files(7).name);
% load('location_key');
% 
% for nn = 1:length(PhysPath.phlo)
%     temp_locs = PhysPath.phlo(nn).all;
%   for aa = 1:length(temp)  
%     for bb = 1:length(location_key)
%         if strcmp([temp_locs{aa, 1}], [location_key{bb,1}])==1
%         tloc_num(aa, 1) = location_key{bb, 2};
%         end
%     end
%   end
% end
% 
% %%
% 
% [B, ~, ib] = unique(location_number, 'rows');
% numoccurrences = accumarray(ib, 1);
% indices = accumarray(ib, find(ib), [], @(rows){rows});  %the find(ib) simply generates (1:size(a,1))'
% %%
% grandcorr = zeros(26, 26);
% for ii = 1:length(B)
%     grandcorr(B(ii, 1), B(ii, 2)) = numoccurrences(ii);
% end
% 
% fig=figure_ctrl('IED-SPI pairs', 1000,1000); imagesc(grandcorr'); colormap jet; axis xy;
% savefig('IEDSPI_pairs.fig');
% savejpg(fig);
%         
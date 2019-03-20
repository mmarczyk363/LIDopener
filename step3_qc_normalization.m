clc; clearvars; close all
addpath('scripts')

% load mat files
load('data\data_match.mat')
n = size(data_norm,2);

% normalize using upper-quantile method
data_raw = [data_norm, data_canc];
[data,fea,~,del] = uq_norm(data_raw,fea,zeros(1,length(fea)));
saveas(gcf,'res\QC\norm_check.fig'); saveas(gcf,'res\QC\norm_check.png'); close gcf;

% save normalized data
data_norm = log2(data(:,1:n)+1);
data_canc = log2(data(:,(n+1):end)+1);
save('data\data_match_norm.mat','data_norm','data_canc','fea','patient')

% calculate PCA
[~,score_raw] = pca(log2(data_raw+1)');
[~,score] = pca(log2(data+1)');

group = false(1,size(data,2));
group((n+1):end) = true;

%Healthy vs Cancer
figure; hold on; box on; 
plot(score_raw(group==0,1),score_raw(group==0,2),'g.',score_raw(group==1,1),score_raw(group==1,2),'r.')
legend({'Normal','Cancer'})
title('Raw data'); xlabel('1st Principal Component'); ylabel('2nd Principal Component');
saveas(gcf,'res\QC\pca2D_raw.fig'); saveas(gcf,'res\QC\pca2D_raw.png'); close gcf;

figure; hold on; box on; 
plot(score(group==0,1),score(group==0,2),'g.',score(group==1,1),score(group==1,2),'r.')
legend({'Normal','Cancer'})
title('Normalized data'); xlabel('1st Principal Component'); ylabel('2nd Principal Component');
saveas(gcf,'res\QC\pca2D_norm.fig'); saveas(gcf,'res\QC\pca2D_norm.png'); close gcf;

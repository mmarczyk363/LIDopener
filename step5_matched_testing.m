clc; clearvars; close all;
addpath('scripts')

% name of the ENZYME database
enz_name = 'enz_5_mult';

% load expression data with enzyme info
load(['data\data_match_norm_' enz_name '.mat'])

sig_lev = 0.05;
enz_nb = length(enzymes_expr.Name);
isoz_nb = size(enzymes_expr_isos,1);
p_isoz = zeros(isoz_nb,1);
eta2_isoz = p_isoz;
l2FC = zeros(isoz_nb,1);
for b=1:isoz_nb
    clc;disp(['Iter no. ' num2str(b) '/' num2str(isoz_nb)])

    %fit repeated measures model
    data_table = table;
    data_table.normal = data_enz_norm(b,:)';
    data_table.cancer = data_enz_canc(b,:)';
    within = table([{'normal'},{'cancer'}]','VariableNames',{'Condition'});
    rm = fitrm(data_table,'cancer-normal ~ 1','WithinDesign',within,'WithinModel','separatemeans');

    %perform repeated measures anova
    tbl = mauchly(rm);
    ranovatbl = ranova(rm);
    p_isoz(b) = ranovatbl.pValueGG(1);
    eta2_isoz(b) = ranovatbl.SumSq(1)/sum(ranovatbl.SumSq);

    %calculate log2 fold change
    l2FC(b)= (mean(data_enz_canc(b,:)) - mean(data_enz_norm(b,:)));
end

%BH FDR correction
p_isoz_corr = mafdr(p_isoz,'BHFDR',true);
disp(['Significant change with condition: ' num2str(sum(p_isoz_corr(:,1)<sig_lev))])

%plot distributions of p-values
figure; subplot(3,1,1);hist(p_isoz,0.025:0.05:0.975); xlabel('P-value');
subplot(3,1,2);hist(p_isoz_corr,0.025:0.05:0.975); xlabel('Corrected p-value');
subplot(3,1,3);plot(log10(p_isoz),log10(p_isoz_corr),'.');xlabel('P-value'); ylabel('Corrected p-value')
saveas(gcf,'res\DEGs\p_val_dist.fig'); saveas(gcf,'res\DEGs\p_val_dist.png'); close gcf;

% create volcano plot
mavolcanoplot(data_enz_norm, data_enz_canc, p_isoz_corr,'Labels', fea_enz, 'PlotOnly', true)
xlabel('Log Fold change')
saveas(gcf,'res\DEGs\volcano.fig'); saveas(gcf,'res\DEGs\volcano.png'); close gcf;

save(['res/DEGs/uni_' enz_name '.mat'],'p_isoz*','eta2_isoz','l2FC');
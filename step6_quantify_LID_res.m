clc; clearvars; close all;
addpath('scripts')

sig_lev = 0.05;
enz_name = 'enz_5_mult';

% load data
load(['data/data_match_norm_' enz_name '.mat'])
load(['res/DEGs/uni_' enz_name '.mat'])

n_enz = length(enzymes_expr.ID);
n_iso = length(enzymes_expr_isos.gene);
enz_manova_all_isos = cell(n_enz,1);
enz_manova_all = false(n_enz,1);

% find enzymes with only one increased/not change isozyme in cancer
isoz_manova_decr = p_isoz_corr<sig_lev & l2FC < 0;
enz_tmp = unique(enzymes_expr_isos.ind(isoz_manova_decr));
for b=1:length(enz_tmp)
   enz_ind = find(enzymes_expr_isos.ind == enz_tmp(b)); 
   enz_manova_all_isos{enz_tmp(b)} = isoz_manova_decr(enz_ind);
   if sum(enz_manova_all_isos{enz_tmp(b)}) == length(enz_ind)-1
       enz_manova_all(enz_tmp(b)) = true;
   end
end
res_manova_isos_chosen = tabulate(enzymes_expr.iso_nb(enz_manova_all));
disp(['For ' num2str(sum(enz_manova_all)) ' enzymes one isozyme has increased or not changed expression in cancer and others have decreased expression.'])

% prioritize target isoenzymes
LID_res = find_targets(enz_manova_all,enz_manova_all_isos,...
    data_enz_norm,data_enz_canc,enzymes_expr,enzymes_expr_isos,fea_enz);
disp(['No. of all target enzymatic reactions: ' num2str(size(LID_res,1))])
disp(['No. of prioritized target enzymatic reactions: ' num2str(sum(LID_res.Cluster < 3))])
saveas(gcf,'res/LID/MDS.fig'); saveas(gcf,'res/LID/MDS.png'); close gcf;
saveas(gcf,'res/LID/Clust_K_sel.fig'); saveas(gcf,'res/LID/Clust_k_sel.png'); close gcf;
save(['res/LID/LID_res_' enz_name '.mat'],'LID_res','enz_manova_all*');
writetable(LID_res,['res/LID/LID_res_' enz_name '.txt'],'Delimiter','\t')

% Plot results for all prioritized targets
sel_enzymes = LID_res(LID_res.Cluster < 3,[1,4]);
plot_enzymes(data_enz_norm,data_enz_canc,enzymes_expr,enzymes_expr_isos,sel_enzymes)

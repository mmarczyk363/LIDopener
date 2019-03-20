clc; clearvars; close all
addpath('scripts')

% load enzyme database
enz_thr = 5;    %max. number of isozymes
rel_date = '10_2018'; %release date
load(['data/enzymes_' num2str(enz_thr) '_' rel_date '_mult.mat'])
enz_name = ['enz_' num2str(enz_thr) '_mult'];

% load normalized and matched patients expression data
load('data\data_match_norm.mat')
    
%check common enzymes in expression data and ENZYME database
fea_comm = intersect(fea,enzymes_filt_isos.gene);
tmp = [];   % find all enzymes with selected/common genes
for b=1:length(fea_comm)
    tmp = [tmp;enzymes_filt_isos.ind(strcmp(fea_comm{b},enzymes_filt_isos.gene))];
end
tmp_un = unique(tmp);
tmp_un_iso_nb = zeros(size(tmp_un));
for b=1:length(tmp_un)
    tmp_un_iso_nb(b) = sum(tmp_un(b)==tmp);
end
res_expr = tabulate(tmp_un_iso_nb)';
disp('Distribution of no. of isozymes per enzyme')
disp(res_expr)

%plot distribution of no. of isozymes after matching
figure; bar(res_expr(2,:))
xlabel('No. of isozymes'); ylabel('No. of enzymes')
xlim([0,6]); ylim([0,max(res_expr(2,:))+0.1*max(res_expr(2,:))])
hold on; plot([1.5,1.5],ylim,'r')
text(res_expr(1,:),res_expr(2,:),num2str(res_expr(2,:)'),'vert','bottom','horiz','center'); 
saveas(gcf,'res\QC\no_iso_match.fig'); saveas(gcf,'res\QC\no_iso_match.png'); close gcf;

%remove enzymes with one gene
tmp_un(tmp_un_iso_nb<2) = [];
tmp_un_iso_nb(tmp_un_iso_nb<2) = [];

%create new enzyme annotation for expression data
enzymes_expr = struct;
names = fieldnames(enzymes_filt);
for b=1:length(names)
    enzymes_expr.(names{b}) = enzymes_filt.(names{b})(tmp_un);
end
for b=1:length(enzymes_expr.ID)
    del = false(1,enzymes_expr.iso_nb(b));
    for c=1:enzymes_expr.iso_nb(b)
        del(c) = ~logical(sum(strcmp(fea_comm,enzymes_expr.gene{b}{c})));
    end
    enzymes_expr.gene{b}(del)= [];
    enzymes_expr.uniprot{b}(del)= [];
    enzymes_expr.swissprot{b}(del)= [];
    enzymes_expr.iso_nb(b) = length(enzymes_expr.gene{b});
end

%plot distribution of no. of isozymes
res_hum = tabulate(enzymes_expr.iso_nb)';
figure; bar(res_hum(2,:))
xlabel('No. of isozymes'); ylabel('No. of enzymes')
xlim([0,6]); ylim([0,max(res_expr(2,:))+0.1*max(res_expr(2,:))])
hold on; plot([1.5,1.5],ylim,'r')
text(res_hum(1,:),res_hum(2,:),num2str(res_hum(2,:)'),'vert','bottom','horiz','center'); 
saveas(gcf,'res\QC\no_iso_final.fig'); saveas(gcf,'res\QC\no_iso_final.png'); close gcf;

%create list of all isozymes for expression data
warning off;
enzymes_expr_isos = table;
data_enz_norm = zeros(sum(enzymes_expr.iso_nb),size(data_norm,2));
data_enz_canc = data_enz_norm;
count = 1;
for b=1:length(enzymes_expr.ID)
    clc; disp(['It. ' num2str(b) '/' num2str(length(enzymes_expr.ID))])
    for c=1:enzymes_expr.iso_nb(b)
        enzymes_expr_isos.uniprot{count,1} = enzymes_expr.uniprot{b}{c};
        enzymes_expr_isos.swissprot{count,1} = enzymes_expr.swissprot{b}{c};
        enzymes_expr_isos.gene{count,1} = enzymes_expr.gene{b}{c};
        enzymes_expr_isos.ID(count,1) = enzymes_expr.ID(b);
        enzymes_expr_isos.ind(count,1) = b;
        ind = find(strcmp(enzymes_expr_isos.gene{count,1},fea));
        data_enz_norm(count,:) = data_norm(ind,:);
        data_enz_canc(count,:) = data_canc(ind,:);
        count = count + 1;
    end
end
warning on;
disp(['No. of unique isozymes: ' num2str(length(unique(enzymes_expr_isos.gene)))])
disp(['No. of enzymes: ' num2str(length(enzymes_expr.ID))])

fea_enz = enzymes_expr_isos.gene;
group_enz = enzymes_expr_isos.ind;

save(['data\data_match_norm_' enz_name '.mat'],'data_enz*','fea_enz','group_enz',...
    'enzymes_expr','enzymes_expr_isos','patient')
clc; clearvars; close all
addpath('scripts')

% load mat file with all data
load('data\data_all.mat')
    
% Match cancer and normal samples for the same patients
pat_comm = intersect(patient(group==1),patient(group==0));
data_norm = nan(length(fea),length(pat_comm));
data_canc = data_norm;
for b=1:length(pat_comm)
    tmp = find(strcmp(pat_comm{b},patient));
    if length(tmp) ~= 2 % more than 2 samples for the same patient
        disp([num2str(length(tmp)) ' samples for ' pat_comm{b}])
        disp(barcode(tmp))
        tmp2 = group(tmp);
        del = false(size(tmp2));
        for c=0:1
            if sum(tmp2==c) > 1
                tmp3 = find(tmp2==c);
                tmp4 = find(tmp2~=c);
                edit_dist = zeros(size(tmp3));  % calculate edit distance between barcodes
                for d=1:length(tmp3)
                    tmp_dist = zeros(size(tmp4));
                    for e=1:length(tmp4)
                         tmp_dist(e) = editdist(barcode{tmp(tmp4(e))},barcode{tmp(tmp3(d))});
                    end
                    edit_dist(d) = mean(tmp_dist);
                end
                [~,ind_dist] = sort(edit_dist); % select the one with the smallest edit distance of the barcode
                del(tmp3(ind_dist(2:end))) = true;
            end
        end
        tmp(del) = [];
        disp('Remained samples:')
        disp(barcode(tmp))
    end
    data_norm(:,b) = data(:,tmp(~group(tmp)));
    data_canc(:,b) = data(:,tmp(group(tmp)));
end

% find and remove genes with no expression after data matching
del = sum([data_norm,data_canc] == 0,2) == 2*size(data_norm,2);
disp([num2str(sum(del)) ' no expression features removed.']);
data_norm(del,:) = [];
data_canc(del,:) = [];
fea(del) = [];

%check low-expression genes in both groups
nan_canc = 100*sum(data_canc == 0,2)/size(data_canc,2);
nan_norm = 100*sum(data_norm == 0,2)/size(data_norm,2);
figure; subplot(2,1,1)
hist_step = 10;
hist(nan_norm(nan_norm>0),hist_step/2:hist_step:100); xlabel('% of samples with no expr. in normal'); ylabel('No. of genes')
subplot(2,1,2)
hist(nan_canc(nan_canc>0),hist_step/2:hist_step:100); xlabel('% of samples with no expr. in cancer'); ylabel('No. of genes')
saveas(gcf,'res\QC\no_expr_dist.fig'); saveas(gcf,'res\QC\no_expr_dist.png'); close gcf;
figure; hist(nan_canc-nan_norm,-70+hist_step:hist_step:70)
xlim([-70-hist_step,70+hist_step])
xlabel('Difference of no expr. samples between cancer and normal'); ylabel('No. of genes')
saveas(gcf,'res\QC\no_expr_diff.fig'); saveas(gcf,'res\QC\no_expr_diff.png'); close gcf;

figure; hold on; box on;
hist_thr = 40; hist_thr2 = 25;
scatter(nan_norm,nan_canc,'.'); 
plot([0,100],[0,100],'r')
plot([hist_thr,100],[0,100-hist_thr],'r--')
plot([0,100-hist_thr],[hist_thr,100],'r--')
rectangle('Position',[hist_thr2,hist_thr2,100-hist_thr2,100-hist_thr2])
xlabel('% of samples with no expr. in normal'); ylabel('% of samples with no expr. in cancer')
set(datacursormode(gcf),'UpdateFcn',@(obj,event_obj)MyCallback(obj,event_obj,fea,[nan_norm,nan_canc]))
saveas(gcf,'res\QC\no_expr_comp.fig'); saveas(gcf,'res\QC\no_expr_comp.png'); close gcf;

%remove genes with > 25% of no expression samples in both conditions, except those with |difference| higher than 40%.
del = nan_canc > hist_thr2 & nan_norm > hist_thr2 & abs(nan_canc-nan_norm) < hist_thr;
data_canc(del,:) = [];
data_norm(del,:) = [];
fea(del) = [];
disp([num2str(sum(del)) ' genes with > ' num2str(hist_thr2) '% of no expression samples in both conditions, except those with |difference| higher than ' num2str(hist_thr) '%.'])
disp([num2str(length(fea)) ' genes left'])

%update names of remaining genes to the same annotation as in ENZYME
%database
load('data/genes_table.mat')
n_changed = 0;
for b=1:length(fea)
    tmp = sum(strcmp(fea{b},fea_new));
    if tmp<1
        disp([fea{b} ' replaced with ' fea_new{strcmp(fea{b},fea_old)}])
        fea{b} = fea_new{strcmp(fea{b},fea_old)};
        n_changed = n_changed + 1;
    end
end
disp([num2str(n_changed) ' gene names updated.'])

patient = pat_comm;
save('data\data_match.mat','data_norm','data_canc','fea','patient')
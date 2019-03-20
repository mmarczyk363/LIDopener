function table_res = find_targets(enz_chosen,enz_chosen_isos,data_enz_norm,data_enz_canc,enzymes_expr,enzymes_expr_isos,fea_enz,r_tmp)

ind = find(enz_chosen);
%calculate stat. diff. between isozymes
p = zeros(length(ind),1);
for a=1:length(ind)
        ind2 = find(enzymes_expr_isos.ind == ind(a));
        data_test = data_enz_canc(ind2,:) - data_enz_norm(ind2,:); 
        if size(data_test,1) == 2
            [~,p(a)] = ttest(data_test(2,:) - data_test(1,:));
        else
            p(a) = anova1(data_test',[],'off');
        end
end

% calculate ratio of iso1/iso_other in normal and cancer
l2FC_iso = zeros(length(ind),3);
target_expr_mean = zeros(length(ind),2);
other_expr_max = target_expr_mean;
target_expr_descr = cell(length(ind),3);
for a=1:length(ind)
    ind_iso = find(enzymes_expr_isos.ind == ind(a));
    target_expr_descr{a,1} = enzymes_expr.ID{ind(a)};
    target_expr_descr{a,2} = enzymes_expr.Name{ind(a)};
    target_expr_descr{a,3} = enzymes_expr.iso_nb(ind(a));
    target_expr_descr{a,4} = fea_enz{ind_iso(~enz_chosen_isos{ind(a)})};
    target_expr_descr{a,5} = enzymes_expr_isos.uniprot(ind_iso(~enz_chosen_isos{ind(a)}));
    if target_expr_descr{a,3} < 3
        target_expr_descr{a,6} = fea_enz{ind_iso(enz_chosen_isos{ind(a)})};
    else
        tmp1 = fea_enz(ind_iso(enz_chosen_isos{ind(a)}));
        target_expr_descr{a,6} = tmp1{1};
        for b=2:target_expr_descr{a,3}-1
            target_expr_descr{a,6} = strcat(target_expr_descr{a,6},{','},tmp1{b});
        end
    end
    target_expr_mean(a,1) = mean(data_enz_norm(ind_iso(~enz_chosen_isos{ind(a)}),:));
    target_expr_mean(a,2) = mean(data_enz_canc(ind_iso(~enz_chosen_isos{ind(a)}),:));
    other_expr_max(a,1) = max(mean(data_enz_norm(ind_iso(enz_chosen_isos{ind(a)}),:),2));
    other_expr_max(a,2) = max(mean(data_enz_canc(ind_iso(enz_chosen_isos{ind(a)}),:),2));
    l2FC_iso(a,1) = target_expr_mean(a,1) - other_expr_max(a,1);
    l2FC_iso(a,2) = target_expr_mean(a,2) - other_expr_max(a,2);
    l2FC_iso(a,3) = l2FC_iso(a,2) - l2FC_iso(a,1);
end

%select the best candidates
score_all = [tiedrank(target_expr_mean(:,2)),tiedrank(-other_expr_max(:,2)),...
    tiedrank(-abs(l2FC_iso(:,1))),tiedrank(l2FC_iso(:,2)),tiedrank(l2FC_iso(:,3))];
score = sum(score_all,2);
score = 100*(score - min(score))/range(score);

%cluster measures
[idx,c_ind] = kmeans_cluster(score_all,20,'sqeuclidean');
idx = idx(:,c_ind(3));  %choose the no. of clusters by gap statistic
cl_score = zeros(1,length(unique(idx)));
for a=1:length(unique(idx))
    cl_score(a) = median(score(idx == a));
end
[~,idx_sort] = sort(cl_score,'descend'); 
idx = idx + 20;
for a=1:length(unique(idx))
    idx(idx==20+idx_sort(a)) = a;
end

%save results in one table
table_res = [cell2table(target_expr_descr,'VariableNames',{'ID','Name','Isos','GeneName','Protein','GeneNameOther'}),...
    array2table([target_expr_mean,other_expr_max,l2FC_iso,p,score,idx,score_all],'VariableNames',{'TargetNormal','TargetCancer',...
    'OtherNormal','OtherCancer','DiffNormal','DiffCancer','Change','PValIso','Score','Cluster','M1','M2','M3','M4','M5'})];
table_res = sortrows(table_res,{'Cluster','Score'},{'ascend','descend'});

%visualize clustering results
dissimilarities = pdist(score_all,'cityblock'); % Create a dissimilarity matrix.
score = mdscale(dissimilarities,2,'Start','random');
sel = idx;
sel(sel>2) = 0;
figure; plot(score(sel==1,1),score(sel==1,2),'r*',score(sel==2,1),score(sel==2,2),'ro',score(sel==0,1),score(sel==0,2),'k.')
legend({'Cluster 1','Cluster 2','NS'})
xlabel('MDS 1'); ylabel('MDS 2')

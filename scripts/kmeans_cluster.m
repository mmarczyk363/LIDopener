function [idx,c_ind,stats] = kmeans_cluster(data,K,method)

%define basic variables
B = 100;
c_ind = zeros(1,3);
rep_no = 100;

%calculate matrix of distances
c_dist = pdist(data);

idx = nan(size(data,1),K); 
stats = struct();
s = nan(1,K); DI = s; W = s; 
W_est = nan(B,K);
for a=2:K
%     disp(['Number of clusters: ' num2str(a) '/' num2str(K)])
    rng(3); % For reproducibility
    idx(:,a) = kmeans(data,a,'Replicates',rep_no,'Distance',method);
    s(a) = mean(silhouette(data,idx(:,a)));
    DI(a) = dunn_index(c_dist,idx(:,a));
    W(a) = w_calc(c_dist,idx(:,a));
    W_est(:,a) = w_est_calc(data,a,B);
end

%find optimal no. of clusters
[~,c_ind(1)] = nanmax(s);   %silhouette method
[~,c_ind(2)] = nanmax(DI);  %dunn index
gap = nan(1,K); sd = gap; sk = gap; %gap statistic
for a=2:K
    gap(a) = sum(log(W_est(:,a))-log(W(a)))/B;
    sd(a) = sqrt(sum(log(W_est(:,a)-sum(log(W_est(:,a)))/B).^2)/B);
    sk(a) = sd(a)*sqrt(1+1/B);
end
tmp = find(diff(gap) < 0);
if ~isempty(tmp); c_ind(3) = tmp(1); else; c_ind(3) = K; end

figure; subplot(3,1,1); hold on; box on;
plot(2:K,s(2:end)); ylabel('Mean silhouette')
plot([c_ind(1),c_ind(1)],ylim,'r')
subplot(3,1,2); hold on; box on;
plot(2:K,DI(2:end)); ylabel('Dunn index')
plot([c_ind(2),c_ind(2)],ylim,'r')
subplot(3,1,3); hold on; box on;
plot(2:K,gap(2:end)); ylabel('Gap statistic')
plot([c_ind(3),c_ind(3)],ylim,'r')
xlabel('No. of clusters'); 

stats.s = s;
stats.W = W;
stats.DI = DI;
stats.gap = gap;
stats.sk = sk;
stats.idx = idx;

function W_est = w_est_calc(data,k,B)

[n,p] = size(data);
W_est = zeros(1,B);
for a=1:B
    data_temp = zeros(n,p);
    for c=1:p
        data_temp(:,c) = min(data(:,c)) + range(data(:,c)).*rand(n,1);
    end
    idx = kmeans(data_temp,k);
    c_dist = pdist(data_temp);
    W_est(a) = w_calc(c_dist,idx);
end

function W = w_calc(c_dist,ind)

c_dist = squareform(c_dist);
clust = unique(ind);
n = length(clust);
W = zeros(1,n);
for a=1:n
    idx = find(ind==clust(a));
    m = length(idx); suma = 0;
    for b =1:m
        for c=1:m
            suma = suma + c_dist(idx(b),idx(c));
        end
    end
    W(a) = suma/(2*m);
end
W = sum(W);

function DI=dunn_index(c_dist,ind)

c_dist = squareform(c_dist);
c_numb = length(unique(ind));
d_inter = zeros(1,c_numb); 
d_intra = d_inter;
DI = d_inter;

for a=1:c_numb
    d_inter(a) = max(max(c_dist(ind==a,ind==a)));
    d_intra(a) = max(max(c_dist(ind==a,ind~=a)));
end

d_max = max(d_inter);
ind2 = 1:c_numb;
for a=1:c_numb
    DI(a) = min(d_intra(ind2~=a)/d_max);
end
DI = min(DI);
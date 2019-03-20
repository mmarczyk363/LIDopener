function [data,fea_names,scale,del] = uq_norm(data,fea_names,del,scale)
%Upper quartile normalization
%
%Implemented by:
%Michal.Marczyk@polsl.pl

%remove features with no expression
n = size(data,2);
if nargin < 3
    del = sum(data==0,2)==n;
%     del = sum(data<6 ,2)>=0.5*n;
end
if sum(del)
    disp([num2str(sum(del)) ' no expression features removed.']);
    data(del,:) = [];
    fea_names(del) = [];
end

%find Q3 and normalize each sample
sum1 = sum(data,1);
sum2 = sum(data>0,1);
q3 = zeros(1,n);
for a=1:n
    q3(a) = quantile(data(:,a),.75);
    if (q3(a) == 0)
        disp('Third quantile equals 0.')
    end
    data(:,a) = data(:,a)/q3(a);
end
figure; subplot(3,1,1);
plot(q3); xlabel('Sample'); ylabel('Q3')
subplot(3,1,2);
plot(sum1,q3,'.'); lsline()
xlabel('Sum of counts'); ylabel('Q3')
subplot(3,1,3);
plot(sum2,q3,'.'); lsline()
xlabel('Sum of expressed genes'); ylabel('Q3')

if nargin<4
    scale = mean(q3);
end
data = data * scale; %scale up data
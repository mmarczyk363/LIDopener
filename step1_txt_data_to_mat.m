clc; clearvars; close all

% load txt data
d = dir('data'); d(1:2) = [];
n = length(d);

% get names of the data (cancer samples set must end with _cancer and normal samples with _normal)
names = cell(n,1);
for a=1:n
    tmp = strsplit(d(a).name,'.txt');
    names{a} = tmp{1};
end
T1 = readtable(['data\',names{1},'.txt']);
T2 = readtable(['data\',names{2},'.txt']);

% match gene names
if sum(strcmp(table2cell(T1(:,1)),table2cell(T2(:,1)))) == size(T1,1)
    disp('Matching ok.')
    
    % merge data into one array
    data = [table2array(T1(:,3:end)),table2array(T2(:,3:end))];
    fea = table2cell(T1(:,1));
    barcode = [T1.Properties.VariableNames(3:end),T2.Properties.VariableNames(3:end)]';
    patient = cell(size(barcode));
    for b=1:length(barcode)
        tmp = find(barcode{b}=='_');
        barcode{b}(tmp) = '-';
        patient{b} = barcode{b}(1:tmp(3)-1);
    end
    

    % annotate normal and cancer samples as group
    group_names = [repmat(names(1),1,size(T1,2)-2),repmat(names(2),1,size(T2,2)-2)];
    group = false(1,size(data,2));
    group(contains(group_names,'cancer')) = true;
    disp(['Matched patients: ' num2str(length(intersect(patient(group==0),patient(group==1))))])
    
    % save dataset
    save('data\data_all.mat','data','fea','group','group_names','barcode','patient');
end
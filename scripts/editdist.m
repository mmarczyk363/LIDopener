function [distance, varargout] = editdist(str1, str2)
%EDITDIST computes the mininum edit distance between two strings
%   this metric is used to determine the dissimilarity between sequences

cost = 1;                                       % delete and insert cost
sub = 1;                                        % substitute cost

% initialize the grid with extra row and column for null strings
D = zeros(length(str1)+1, length(str2)+1);      % add extra row and column
D(:,1) = 0:length(str1);                        % the null string
D(1,:) = 0:length(str2);                        % the null string

% loop through a pair of characters
for i = 2:size(D, 1)                            % skip the null row
    for j = 2:size(D, 2)                        % skip the null col
        deletion = D(i-1, j) + cost;            % compute deletion cost
        insertion = D(i, j-1) + cost;           % compute insertion cost
        if str1(i-1) == str2(j-1)               % if the same character
            substitution = D(i-1, j-1) + 0;     % no substitution cost 
        else
            substitution = D(i-1, j-1) + sub;   % compute substitution cost
        end
        D(i,j) = min([deletion, insertion,...   % use the minimum of three
            substitution]);    
    end
end
% finish the computation
distance = D(end);                              % return the final element
rownames = regexp(str1, '\w', 'match');         % parse string into chars
rownames = [{'Null'}, lower(rownames)];         % add null 
colnames = regexp(str2, '\w', 'match');         % parse string into chars
colnames = [{'Null'}, lower(colnames)];         % add null 
varargout = {D, rownames, colnames};            % optional return values

end


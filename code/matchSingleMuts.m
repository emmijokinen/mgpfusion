function match=matchSingleMuts(mutE,mutR)
% Find indices of single mutations from mutR that correspond to mutations
% in mutE. Zero, if multiple mutation in mutE or no corresponding mutation
% in mutR.

% INPUT mutE = cell array of mutations (size NE x 1)
%       mutR = cell array of mutations
%
% OUTPUT match = list of indices (size NE x 1)
%

singleE=find(cellfun(@(mut) sum(isletter(mut))==2,mutE));
mutMatching=[singleE,zeros(length(singleE),1)];
j=1;
for i=1:length(singleE)
    mut=mutE{singleE(i)};
    ji=find(cellfun(@(x) strcmp(x,mut),mutR(j:end)),1,'first')+j-1;
    if ~isempty(ji)
        j=ji;
        mutMatching(i,2)=j;
    end
end
match=zeros(size(mutE));
match(mutMatching(:,1))=mutMatching(:,2);
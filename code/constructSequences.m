function sequencesNew=constructSequences(sequenceWT, mutations1, includeWT)
%% CONSTRUCTSEQUENCES creates a matrix of sequences 
% based on the given WT sequence and mutations
%
% INPUT: sequenceWT = sequence of the WT protein
%        mutations = cell array of mutations. One cell contains all
%           mutations that are done to one protein.
%
% OUTPUT: sequences = matrix with one sequence on every row. The first
%           one is the WT protein which is followed by the mutated proteins
%           in the same order as given in mutations
% 

if isletter(sequenceWT(1))
    sequenceWT=double(aa2int(sequenceWT));
end

numMut=numel(mutations1);
if includeWT; wt=1; else; wt=0; end
sequencesNew=repmat(sequenceWT,numMut+wt,1);

for i=1:numMut
    mut=mutations1{i};
    letter_locs=isletter(mut);
    if mod(sum(letter_locs),2)~=0 || any(~ismember(mut,'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'))
        error('Invalid mutation number %d: %s',i,mut)
    end
    letter_inds=find(letter_locs);
    for j=1:2:sum(letter_locs)
        num=str2double(mut(letter_inds(j)+1:letter_inds(j+1)-1));
        if isnan(num) 
            warning('Mutation %s at index %d is not known. The WT sequence is left at this location.',mut,i);
            break
        elseif sequenceWT(num)~=aa2int(mut(letter_inds(j)))
            warning('Mutation %s at index %d does not match the WT, the original AA is %c',mut,i, int2aa(sequenceWT(num)));
        end
        sequencesNew(i+wt,num)=double(aa2int(mut(letter_inds(j+1))));
    end
end
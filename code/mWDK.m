function K=mWDK(WTSeq,mutations,al,depth,normalize,S)
%% MWDK calculates a graph kernel matrix with fixed adjacency matrix
%
% INPUT:    WTSeq = sequence of the wild type protein. (chars or integers)
%           mutations = cell array of strings containing the mutations in
%               format AiB, where A is the amino acid in the WTSeq and B is
%               the mutated residue. i is the residue number. Multiple
%               mutations where there are several of these triplets are
%               also possible (eg. M1AL5I).
%           al = adjacency list. Cell array that contains an array of
%               neighbours for each residue.
%           depth = depth of the neighbourhood. If 1, then only residues
%               in contact are considered, i
%           normalize = boolean. If true, the matrix is normalized
%           S = subsitution matrix (PSD)
%
% OUTPUT:   K = the calculated kernel matrix

sequences=constructSequences(WTSeq,mutations,1); % includes a sequence for the WT

N=size(sequences,1); % num of sequences
n=length(WTSeq); % num of AAs
numNbs=cellfun(@length,al);

% Determine the neighbourhood residues
if depth>1
    alit=cell(n,1);
    for d=1:depth
        for i=1:n
            ali=al{i};
            tmp=al(ali);
            alit{i}=unique([tmp{:,:}]);
        end
        al=alit;
    end
end

K = zeros(N,N);
neighpos = find(numNbs > 0); % no lonely residues
for i=1:length(neighpos)
    r=neighpos(i);
    zs = sequences(:,r); % selectors
    ns = sequences(:,al{r});
    zscores = S(zs,zs);
    nscores = squeeze( sum( bsxfun(@(i,j) S(i,j), ns, permute(ns, [3 2 1])), 2) );
    K = K + zscores .* nscores;
end

if normalize
    K=K./sqrt(diag(K)*diag(K)'); % normalize
end


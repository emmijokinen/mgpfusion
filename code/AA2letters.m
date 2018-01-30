function letters=AA2letters(residues,replace)
%% AA2LETTERS transforms a cell array of strings to a cell array of letters.
% unknown strings (such as ligands) are kept as they were
%
% INPUT:    residues = three letter presentation of amino acids
% OUTPUT:   letters = one letter presentation of amino acids (unknown names 
%               are kept as they were

% MSE is considered as MET


if replace
    % MSE-MET, TRN-TRP, PCA-GLU,CSO-CYS,M3L-LYS
    names=['ALA';'ARG';'ASN';'ASP';'CYS';'GLN';'GLU';'GLY';'HIS';'ILE'; ...
           'LEU';'LYS';'MET';'PHE';'PRO';'SER';'THR';'TRP';'TYR';'VAL';...
           'MSE';'TRN';'PCA';'CSO';'M3L'];
    shorts=['A';'R';'N';'D';'C';'Q';'E';'G';'H';'I';...
            'L';'K';'M';'F';'P';'S';'T';'W';'Y';'V';...
            'M';'W';'E';'C';'K'];
else
    names=['ALA';'ARG';'ASN';'ASP';'CYS';'GLN';'GLU';'GLY';'HIS';'ILE'; ...
           'LEU';'LYS';'MET';'PHE';'PRO';'SER';'THR';'TRP';'TYR';'VAL'];
    shorts=['A';'R';'N';'D';'C';'Q';'E';'G';'H';'I';...
           'L';'K';'M';'F';'P';'S';'T';'W';'Y';'V'];
end

letters=residues;
for i=1:numel(shorts)
    letters(strcmp(names(i,:),residues))={shorts(i)};
end
    
function AA=letters2AA(residues)
%% AA2LETTERS transforms a cell array of strings to a cell array of letters.
% unknown strings (such as ligands) are kept as they were
%
% INPUT:    residues = three letter presentation of amino acids
% OUTPUT:   letters = one letter presentation of amino acids (unknown names 
%               are kept as they were

names=['ALA';'ARG';'ASN';'ASP';'CYS';'GLN';'GLU';'GLY';'HIS';'ILE'; ...
    'LEU';'LYS';'MET';'PHE';'PRO';'SER';'THR';'TRP';'TYR';'VAL';'MSE'];
shorts=['A';'R';'N';'D';'C';'Q';'E';'G';'H';'I';...
    'L';'K';'M';'F';'P';'S';'T';'W';'Y';'V';'M'];

AA=residues;
for i=1:numel(shorts)
    AA(strcmp(shorts(i,:),residues))={names(i,:)};
end
    
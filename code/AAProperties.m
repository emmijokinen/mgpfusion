function structure=AAProperties(structure)
%% AAPROPERTIES returns size, charge and hydrophobicity of given amino acids
% INPUT:    structure that contains fields:
%           sequence.residues = amino acids of the protein
%           numResidues = number of residues in the protein
%           OR
%           structure = sequence of amino acids presented with either one
%           letter or three letter codes
% OUTPUT:   structure = the original structure now also containing fields
%               structure.sequence.sizes, structure.sequence.charge and 
%               structure.sequence.hydropathy
%           OR
%           structure = a numResidues x 3 matrix containing size, charge
%              and hydrophobicity for each residue in the given sequence

% normal amino acids and their properties
sizes=[67;148;96;91;86;114;109;48;118;124;124;135;124;135;90;73;93;163;141;105];
hydropathy=[-1,1,1,1,0, 1,1,-1,0,-1, -1,1,0,-1,-1, 1,1,-1,-1,-1];
charge=[0,1,0,-1,0, 0,-1,0,1,0, 0,1,0,0,0, 0,0,0,0,0];

res={'ALA';'ARG';'ASN';'ASP';'CYS';'GLN';'GLU';'GLY';'HIS';'ILE';...
    'LEU';'LYS';'MET';'PHE';'PRO';'SER';'THR';'TRP';'TYR';'VAL'};

if isstruct(structure)
    if isfield(structure,'numResidues')&& isfield(structure,'sequence')
        numResidues=structure.numResidues;
        residueSeq={structure.sequence.residues};
    else
        error('The structure should have fields numResidues and sequence');
    end
    
    for i=1:numResidues
        ind=find(strcmp(res,residueSeq(i)));
        if isempty(ind)
            if strcmp(residueSeq(i),'MSE')
                % Consider MSE as MET
                structure.sequence(i).size=sizes(13);
                structure.sequence(i).charge=charge(13);
                structure.sequence(i).hydropathy=hydropathy(13);
            else
                % Arbitrary values for ligands (and other non-defined residues
                % which might exist)
                structure.sequence(i).size=200;
                structure.sequence(i).charge=0;
                structure.sequence(i).hydropathy=0;
            end
        else
            structure.sequence(i).size=sizes(ind);
            structure.sequence(i).charge=charge(ind);
            structure.sequence(i).hydropathy=hydropathy(ind);
        end
    end
else
    residueSeq=structure;
    if length(residueSeq{1})==1
        residueSeq=letters2AA(residueSeq);
    end
    
    numResidues=length(structure);
    structure=zeros(numResidues,3);

    for i=1:numResidues
        ind=find(strcmp(res,residueSeq(i)));
        if isempty(ind)
            if strcmp(residueSeq(i),'MSE')
                % Consider MSE as MET
                structure(i,1)=sizes(13);
                structure(i,2)=charge(13);
                structure(i,3)=hydropathy(13);
            else
                % Arbitrary values for ligands (and other non-defined residues
                % which might exist)
                structure(i,1)=200;
                structure(i,2)=0;
                structure(i,3)=0;
            end
        else
            structure(i,1)=sizes(ind);
            structure(i,2)=charge(ind);
            structure(i,3)=hydropathy(ind);
        end
    end
end

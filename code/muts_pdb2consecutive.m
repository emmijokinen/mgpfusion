function muts=muts_pdb2consecutive(mutsPDB,prot)
% changes mutation numbering to that used by Rosetta
% 
% INPUT:    mutsPDB = mutations numbered based on the corresponding PDB
%               structure. Cell of size n x 1
%           prot = structure with protein information
% OUTPUT:   muts = mutations numbered by their position in the amino acid 
%               sequence. Cell of size n x 1
%
    origSeqNum={prot{1}.sequence(:).origSerNo}';
    origSeqNum(cellfun(@isempty, origSeqNum))={nan};
    origSeqNum=cell2mat(origSeqNum);
    
    n=length(mutsPDB);
    properNum=cell(n,1);
    
    nums_pdb=cellfun(@(x) strsplit(x,'[A-Z]+','DelimiterType','RegularExpression'),mutsPDB,'un',0);
    nums_pdb=cellfun(@(x) str2double(x(2:end-1)),nums_pdb,'un',0);
    
    for j=1:n
        properNum(j)=cellfun(@(x) find(ismember(origSeqNum,x)),nums_pdb(j),'un',0);
    end
        
    letter=cellfun(@(x) strsplit(x,'[0-9.]+','DelimiterType','RegularExpression'),mutsPDB,'un',0);
    
    muts=cell(n,1);
    for j=1:n
        try
            tmp=[letter{j};[cellfun(@num2str,num2cell(properNum{j}),'un',0)',{''}]];
            muts{j}=[tmp{:}];
        catch
           disp([mutsPDB{j},' removed. Residue not present in the PDB-structure']) 
        end
    end
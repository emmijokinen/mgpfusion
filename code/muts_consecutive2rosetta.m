function mutsR=muts_consecutive2rosetta(muts,prot)
% changes mutation numbering to that used by Rosetta
% 
% INPUT:    muts = mutations numbered by the position in the amino acid
%               sequence (regardless if they're present in PDB or not).
%               Cell of size n x 1
%           prot = structure with protein information
% OUTPUT:   mutsR = mutations in Rosetta format. Cell of size n x 1
%

    origSerNo={prot{1}.sequence(:).origSerNo};
    I_tmp=(cellfun(@isempty, origSerNo));
    if sum(I_tmp)==0
        mutsR=muts;
    else
        I_cumsum=cumsum(I_tmp);
        I_noempty=~cellfun(@isempty,muts);
        nums_pro=cell(size(muts));
        nums_pro(I_noempty)=cellfun(@(x) strsplit(x,'[A-Z]+','DelimiterType','RegularExpression'),muts(I_noempty),'un',0);
        nums_pro(I_noempty)=cellfun(@(x) str2double(x(2:end-1)),nums_pro(I_noempty),'un',0);
        nums_ros=cellfun(@(x) x-I_cumsum(x),nums_pro,'un',0);
       
        letters=cell(size(muts));
        letters(I_noempty)=cellfun(@(x) strsplit(x,'[0-9]+','DelimiterType','RegularExpression'),muts(I_noempty),'un',0);
    
        mutsR=cell(size(muts));
        for j=1:length(muts)
            try
                tmp=[letters{j};[cellfun(@num2str,num2cell(nums_ros{j}),'un',0),{''}]];
                mutsR{j}=[tmp{:}];
            catch
               disp('error') 
            end
        end         
    end
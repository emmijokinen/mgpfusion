function structure=protein(filename,threshold,useOnlyOneChain,replacements,modelnum)
%PROTEIN processes PDB-structures to get distances, neighbours, and other
%information regarding the residues
%
%   INPUT:  filename = PDB ID from the PDB database or a filename 
%                containing a PDB structure
%           threshold = maximum distance between contacting residues (Å)
%           useOnlyOneChain = set to 0 if all chains should be used.
%              Otherwise set the number or letter of the chain that should 
%              be considered.
%           replacements = determins residue names of (not natural)
%               residues that have been used to replace residues of WT 
%               protein and their corresponding amino acids in the WT
%               protein, give {[]} if none..
%           modelnum = number of the model to be used from the
%               pdb-structure. Default is 1.
%
%   OUTPUT: structure with following fields:
%           contacts = contacts within (and between) the chains
%           numResidues = number of residues in the protein
%           chains.starts = starting indexes of each chain
%           chains.ends = ending indexes of each chain
%           distances = closest distances between residues. If there are no 
%               coordinates, value zero is assigned.
%           sequence = contains information of each residue in the sequence
%           all = array of some information from sequence (X, Y, Z (of 
%               centers), conservations, size, charge, hydropathy, SS 
%               (0: undetermined, 1: Helix, 2: Sheet, 3: Ligand)
%           origSequence = Sequence, where hetatoms are replaced by the 
%               amino acids that occur in the WT sequence
%           score = fraction of Uniprot sequence in pdb
%
%% read pdb and other initalisation
try
    if exist(filename,'file')
        pdbstruct = pdbread(filename);
    else
        pdbstruct = getpdb(filename);
    end
catch
    error('given pdb was not a file or valid pdbid.');
end

% Select model, 1 is default
if nargin<5||isempty(modelnum)
    modelnum=1;
end

% Get sequence information
if isfield(pdbstruct, 'Sequence')
    numSeqs = numel(pdbstruct.Sequence);
else
    numSeqs=numel(unique([pdbstruct.Model(modelnum).Atom.chainID]));
    if numSeqs<1
        numSeqs=1;
    end
end

if isletter(useOnlyOneChain) || useOnlyOneChain
    numSeqs=1;
end

modelstruct = pdbstruct.Model(modelnum);
modelstruct.coords=struct;
modelstruct.order=struct; % order for each sequence
modelstruct.contacts=struct; % contacts for each sequence
modelstruct.contacts2=struct; % contacts between sequences
modelstruct.res=struct; %residues to be given as output

%% Check for missing residues
if isfield(pdbstruct, 'Remark465')
    missing_cell=cellstr(pdbstruct.Remark465);
    for i=1:numel(missing_cell)
        if strfind(missing_cell{i},'RES C SSSEQ')%'M RES C SSSEQ'
            break;
        end
    end
    
    missing=strsplit([missing_cell{i+1:end}]);
    while isempty(missing{1})
        missing(1)=[];
    end
    missing_residues=reshape(missing,3,[])';
    missing_residues(:,3)=num2cell(str2double(missing_residues(:,3)));
else
    missing_residues=[];
end

%% Check DBReferences for UniProt sequence
uSeq=cell(length(pdbstruct.DBReferences),2);
if isfield(pdbstruct, 'DBReferences')
    try
        for i=1:length(pdbstruct.DBReferences)
            if strfind(pdbstruct.DBReferences(i).database,'UNP')
                accessionCode=strtrim(pdbstruct.DBReferences(i).dbAccession);
                fastaR=fastaread(urlread(['http://www.uniprot.org/uniprot/?query=id:' accessionCode '&format=fasta']));
                uSeq{i,1}=pdbstruct.DBReferences(i).chainID;
                uSeq{i,2}=fastaR(1).Sequence;
                score=1;
            end
        end
    catch
        disp('Error downloading URL for checking Uniprot sequence.')
    end
end

if isempty(uSeq)
    score=nan;
    disp([filename ': No uniprot reference / bad reference. Could not determine if the structure it partial']);
end

%% Go through each chain separately
for seq = 1:numSeqs
    % If only one sequence is used, check which one. Otherwise go trough
    % all sequences in order
    if numSeqs==1
        if ischar(useOnlyOneChain)
            seqT=find([pdbstruct.Sequence(:).ChainID]==useOnlyOneChain);
            if ~isempty(seqT)
                seq=seqT; %#ok<FXSET>
            else
                error('No sequence %c for %s',useOnlyOneChain,filename)
            end
        else
            seq=useOnlyOneChain; %#ok<FXSET>
        end
    end
 
    if isfield(pdbstruct, 'Sequence')
        modelstruct.numRes(seq).nr=pdbstruct.Sequence(seq).NumOfResidues;
        % Get the residue names
        modelstruct.res(seq).r=strsplit(pdbstruct.Sequence(seq).ResidueNames)';
        % get IDs for current chain
        chainID = pdbstruct.Sequence(seq).ChainID;
        Ichain = [modelstruct.Atom.chainID]' == chainID;
    else % incase of e.g. swiss model or preminimised structure
        uniqResSeq=unique([modelstruct.Atom.resSeq]);
        modelstruct.numRes(seq).nr=length(uniqResSeq);
        temp=arrayfun(@(n) find(uniqResSeq(n)==[modelstruct.Atom.resSeq],1), 1:modelstruct.numRes(seq).nr);
        modelstruct.res(seq).r=cellstr(reshape([modelstruct.Atom(temp).resName],3,[])');
        chainID=unique([pdbstruct.Model(modelnum).Atom.chainID]);

        if isempty(chainID)
            Ichain= true(length([modelstruct.Atom.resSeq]),1);
            chainID=char('A'-1+seq);
        else
            chainID=chainID(seq);
            Ichain = [modelstruct.Atom.chainID]'==chainID;
        end
    end
    
    % If letters have been useded in the numbering of residues, change them
    % to numbers so that correct order remains: E.g. 27-> 27, 27A -> 27.01, 
    % 27B -> 27.02 and so on, or 1X -> 0.241, 1 -> 1
    needOrderFix=false;
    if isfield(modelstruct.Atom, 'iCode')
        hasIcode=~cellfun(@isempty,{modelstruct.Atom(:).iCode});
        chainHasIcode=hasIcode'&Ichain; 
        if any(chainHasIcode)
            startChain=find(Ichain,1,'first');
            startIcode=find(chainHasIcode,1,'first');
            % First notation, where residue number with letter code comes
            % before the plain number. For these cases the third decimal is one.
            % Second notation (else), is where the plain number comes first
            % and added letters after that. For these cases the third decimal is zero.
            if hasIcode(startChain) || modelstruct.Atom(startIcode-1).resSeq<modelstruct.Atom(startIcode).resSeq
                tmp=num2cell([modelstruct.Atom(chainHasIcode).resSeq]+...
                    (double([modelstruct.Atom(chainHasIcode).iCode])-double('A')+1)/100-0.999);
                [modelstruct.Atom(chainHasIcode).resSeq]=tmp{:};
                needOrderFix=sum(chainHasIcode);           
            else
                tmp=num2cell([modelstruct.Atom(chainHasIcode).resSeq]+...
                    (double([modelstruct.Atom(chainHasIcode).iCode])-double('A')+1)/100);
                [modelstruct.Atom(chainHasIcode).resSeq]=tmp{:};
                needOrderFix=sum(chainHasIcode);
            end
        end
    end
     
    % get secondary structure starts and ends
    modelstruct.secondary(seq).s=cell(modelstruct.numRes(seq).nr,1);
    if isfield(pdbstruct,'Sheet')
        sheetLocations=[pdbstruct.Sheet.initChainID]'==chainID;
        sheets=[[pdbstruct.Sheet(sheetLocations).initSeqNum]' [pdbstruct.Sheet(sheetLocations).endSeqNum]'];
    end
    if isfield(pdbstruct,'Helix')
        helixLocations=[pdbstruct.Helix.initChainID]'==chainID;
        helixes=[[pdbstruct.Helix(helixLocations).initSeqNum]' [pdbstruct.Helix(helixLocations).endSeqNum]'];
    end

    % Find all alpha carbons
    ca=arrayfun(@(n) strcmp(modelstruct.Atom(n).AtomName, 'CA'), 1:numel(modelstruct.Atom));
    
    % If there are heterogen atoms, that are included in the sequence,
    % include them (eg. MSE)
    diff=setdiff(modelstruct.res(seq).r,{modelstruct.Atom(Ichain).resName});
    if ~isempty(diff) && isfield(modelstruct, 'HeterogenAtom') ...
            && ~all(ismember(diff,{'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}))
        hetCoords = ([modelstruct.HeterogenAtom.chainID] == chainID)';
        for i=1:numel(diff)
            hetCoords=hetCoords+~cellfun('isempty',strfind({modelstruct.HeterogenAtom.resName}',diff{i}))>1;
        end

        ca_het = arrayfun(@(n) strcmp(modelstruct.HeterogenAtom(n).AtomName,'CA'),1:numel(modelstruct.HeterogenAtom));
        ca_het = ca_het+hetCoords'>1;
        ca=[ca(Ichain) ca_het(hetCoords)]';
        
        % Create Nx3 matrix of 3D coords        
        modelstruct.coords(seq).c =...
            [modelstruct.Atom(Ichain).X modelstruct.HeterogenAtom(hetCoords).X;
            modelstruct.Atom(Ichain).Y modelstruct.HeterogenAtom(hetCoords).Y;
            modelstruct.Atom(Ichain).Z modelstruct.HeterogenAtom(hetCoords).Z; ]';
        
        % Residues that the the atoms are part of
        modelstruct.order(seq).o = [modelstruct.Atom(Ichain).resSeq,...
            modelstruct.HeterogenAtom(hetCoords).resSeq;];

        % reorder the matrix so that any Het atoms that should be intersperced
        % with normal atoms are in the correct order.
        [modelstruct.order(seq).o, perm] = sort(modelstruct.order(seq).o);
        modelstruct.coords(seq).c = modelstruct.coords(seq).c(perm,:);
        ca=ca(perm,:);
        fprintf('%s: The following heterogen Atoms where added: ',filename);
        disp(diff)
        
        % resSeq: Residue number sequence 
        % strucSeq: residue name sequence of the residues present in the 3D-structure
        ResSeqA=[modelstruct.Atom(Ichain).resSeq];
        uniqResSeqA=unique(ResSeqA);
        resCoordsA=false(size(ResSeqA));
        for i=1:length(uniqResSeqA)
            resCoordsA(find(ResSeqA==uniqResSeqA(i),1,'first'))=true;
        end
        ResSeqH=[modelstruct.HeterogenAtom(hetCoords).resSeq];
        uniqResSeqH=unique(ResSeqH);
        resCoordsH=false(size(ResSeqH));
        for i=1:length(uniqResSeqH)
            resCoordsH(find(ResSeqH==uniqResSeqH(i),1,'first'))=true;
        end
        resCoords=[resCoordsA,resCoordsH];
        
        strucSeq={modelstruct.Atom(Ichain).resName,...
            modelstruct.HeterogenAtom(hetCoords).resName;};
        strucSeq=strucSeq(perm);
        strucSeq=strucSeq(resCoords(perm));
            
        resSeq=[ResSeqA,ResSeqH];
        resSeq=resSeq(perm);  
    else
        % Atom coordinates
        modelstruct.coords(seq).c = [modelstruct.Atom(Ichain).X;
            modelstruct.Atom(Ichain).Y;
            modelstruct.Atom(Ichain).Z; ]';
        % Residues that the atoms are part of
        modelstruct.order(seq).o = [modelstruct.Atom(Ichain).resSeq];
        % Find alpha carbons in the current sequence
        ca=ca(Ichain)';
        
        % resSeq: Residue number sequence 
        % strucSeq: residue name sequence of the residues present in the 3D-structure
        resSeq=[modelstruct.Atom(Ichain).resSeq];
        uniqResSeq=unique(resSeq);
        resCoords=false(size(resSeq));
        for i=1:length(uniqResSeq)
            resCoords(find(resSeq==uniqResSeq(i),1,'first'))=true;
        end
        strucSeq={modelstruct.Atom(Ichain).resName};
        strucSeq=strucSeq(resCoords);
    end
    
    if isfield(pdbstruct, 'Sequence')
        % origSeq: Sequence, where hetatoms are replaced by the amino acids
        % that occur in the WT sequence
        origSeq=(AA2letters(strsplit(pdbstruct.Sequence(seq).ResidueNames),true))';
        tmp=find(cellfun(@(x) length(x)>1,origSeq));
        num_tmp=numel(tmp);
        % if there are residues that were not determined by AA2letters...
        if ~isempty(tmp)
            % if replacements have been defined in the function call..
            if ~isempty(replacements{1,1})
                for i=1:length(tmp)
                    i_tmp=strfind(origSeq{tmp(i)},replacements(:,1));
                    if ~isempty(i_tmp)
                        origSeq(tmp(i))=replacements(i_tmp,2);
                        num_tmp=num_tmp-1;
                        tmp(i)=nan;
                    end
                end
            end
            % otherwise check from UniProt sequence if it was successfully
            % fetched
            if num_tmp>0
                tmp(isnan(tmp))=[];
                old_orig_tmp=origSeq{tmp};
                origSeq(tmp)={'X'};
                origSeq=char(origSeq)';
                if ~isempty(uSeq)
                    uuSeq=uSeq{strcmp(uSeq(:,1),chainID),2};
                    [~,alig]=nwalign(uuSeq,origSeq, 'ExtendGap',0);
                    x_loc=alig(3,:)=='X';
                    origSeq(tmp)=alig(1,x_loc);
                    disp(['Non standard residue ' old_orig_tmp ' was replaced with ' origSeq(tmp)])
                else
                    error(['Non standard residue: ' origSeq{tmp} ' in the sequence. Could not determine origin.'])
                end
            end
        else
            origSeq=char(origSeq)';
        end
        
        % Check if the model is partial, and if yes, how much does it cover
        % from from the complete protein chain (obtained from UniProt)
        if ~isempty(uSeq)
            uuSeq=uSeq{strcmp(uSeq(:,1),chainID),2};
            [~,alig]=nwalign(uuSeq,origSeq, 'ExtendGap',0);
            score=sum(alig(2,:)=='|')/size(alig,2);
            if score<1
                disp([filename ': The structure might be only partial. It contains only ' num2str(score*100) '% of the corresponding Uniprot sequence']);
            end
        end
        
        % Convert strucSeq to letters, replace hetatoms by the amino acids 
        % that occur in the WT sequence
        strucSeq=(AA2letters(strucSeq,true))';
        tmp=find(cellfun(@(x) length(x)>1,strucSeq));
        if ~isempty(tmp)
            if ~isempty(replacements{1,1})
                for i=1:length(tmp)
                    i_tmp=strfind(strucSeq{tmp(i)},replacements(:,1));
                    strucSeq(tmp(i))=replacements(i_tmp,2);
                end
            end
            if num_tmp>0
                tmp(isnan(tmp))=[];
                strucSeq(tmp)={'X'};
                strucSeq=char(strucSeq)';
               
                [~,alig]=nwalign(origSeq,strucSeq, 'ExtendGap',0);
                x_loc=alig(3,:)=='X';
                strucSeq(tmp)=alig(1,x_loc);
            end
        else
            strucSeq=char(strucSeq)';
        end

        % Correct the numbering of residues to go from 1 to number of residues
        [~,alig]=nwalign(strucSeq,origSeq, 'ExtendGap',0);
        % if the numbering is changed, change also secondary structure numbering
        i_s=find(alig(1,:)~='-',1);
        startindex=0;  
        if ceil(modelstruct.order(seq).o(1))~=i_s
            startindex=i_s-ceil(modelstruct.order(seq).o(1));
            modelstruct.order(seq).o=modelstruct.order(seq).o+startindex;
        end 
        strucSeq=alig(1,:);
        
        if needOrderFix
            nonIntegers=mod([modelstruct.order(seq).o],1)~=0;
            uNonIntegers=unique(modelstruct.order(seq).o(nonIntegers));
            for i=1:numel(uNonIntegers)
                select=[modelstruct.order(seq).o]==uNonIntegers(i);
                iFirst=find(select,1,'first');
                if iFirst==1
                    [modelstruct.order(seq).o(select)]=1;
                else
                    [modelstruct.order(seq).o(select)]=modelstruct.order(seq).o(iFirst-1)+1;
                end
                modelstruct.order(seq).o(iFirst:end)=modelstruct.order(seq).o(iFirst:end)+double(~nonIntegers(iFirst:end));
            end      
        end
        
        if isfield(pdbstruct,'Sheet')
            sheets=sheets+startindex;
        end
        if isfield(pdbstruct,'Helix')
            helixes=helixes+startindex;
        end

        % Remove unnecessary gaps
        if numel(strucSeq)<(max([modelstruct.Atom(Ichain).resSeq])+startindex)+numel(strucSeq)-find(isletter(strucSeq),1,'last')+needOrderFix
            m=max([modelstruct.Atom(Ichain).resSeq])+startindex;
            order_unique=unique(modelstruct.order(seq).o(1:sum([modelstruct.Atom.chainID]==chainID)));
            previous_atoms=sum([modelstruct.Atom.chainID]<chainID);
            % Check if there are gaps in the structural sequence
            if (max(order_unique)-min(order_unique)+1)>numel(order_unique)
                i=1;
                while i<numel(order_unique)-1 && order_unique(i)<=(m-1)      
                    % Locate the gaps
                    if order_unique(i+1)-order_unique(i)~=1
                        orig_seq=modelstruct.Atom(previous_atoms+find(modelstruct.order(seq).o==order_unique(i),1)).resSeq;
                        for j=1+orig_seq:(orig_seq+order_unique(i+1)-order_unique(i)-1) 
                            % If the gap is not due to a known missing residue, remove it
                            if isempty(missing_residues)||~ismember(j,[missing_residues{[missing_residues{:,2}]==chainID,3}])
                                 index=find(modelstruct.order(seq).o>order_unique(i),1);
                                 modelstruct.order(seq).o(index:end)=modelstruct.order(seq).o(index:end)-1;
                                 if isfield(pdbstruct,'Sheet')
                                    sheets=sheets+(sheets>=order_unique(i));
                                 end
                                 if isfield(pdbstruct,'Helix')
                                    helixes=helixes+(helixes>=order_unique(i));
                                 end
                                 order_unique=unique(modelstruct.order(seq).o);
                            end
                        end
                    end
                    i=i+1;
                end
            end
        end
    else
        % If there is no sequence field in the pdbstruct, no other adjustments
        % can be done, but to check that the indexing starts with one.
        if modelstruct.order(seq).o(1)~=1
            modelstruct.order(seq).o=modelstruct.order(seq).o-(modelstruct.order(seq).o(1)-1);
        end
        fprintf('%s: The pdb-structure does not contain field Sequence,starting index is set to 1\n',filename);
        origSeq=[];
    end
    modelstruct.score(seq).s=score;
    
    % If the secondary structures are given in the pdb-file, save them
    if isfield(pdbstruct,'Sheet')
        sheetCoords=zeros(modelstruct.numRes(seq).nr,1);
        for i=1:size(sheets,1)
            sheetCoords=sheetCoords+arrayfun(@(n) any(n==sheets(i,1):sheets(i,2)), 1:modelstruct.numRes(seq).nr)';
        end
        modelstruct.secondary(seq).s(logical(sheetCoords))= {'S'};
    end
    if isfield(pdbstruct,'Helix')
        helixCoords=zeros(modelstruct.numRes(seq).nr,1);
        for i=1:size(helixes,1)
            helixCoords=helixCoords+arrayfun(@(n) ...
                any(n==helixes(i,1):helixes(i,2)), 1:modelstruct.numRes(seq).nr)';
        end 
        modelstruct.secondary(seq).s(logical(helixCoords))= {'H'};
    end
    
    % If a residue has several alpha carbons (this should be due to imaging
    % problems) keep only the first one.
    i=1;
    while i<numel(ca)-1
        if ca(i)
            j=i+1;
            while j<numel(ca) && ~ca(j)
                j=j+1;
            end
            if modelstruct.order(seq).o(i)==modelstruct.order(seq).o(j)
               ca(j)=0;
            end
            i=j+1;
        end
        i=i+1;
    end
    
    % Determine the "center" of each residue by the location of its alpha
    % carbon. If alpha carbon is missing, take the mean of other atoms.
    modelstruct.centers(seq).c=cell(1,modelstruct.numRes(seq).nr);
    % Check if ca is missing from residue with other atoms present
    if sum(ca)<length(strucSeq)
        diff=setdiff(resSeq,resSeq(ca));
        for i=1:length(diff)
            modelstruct.centers(seq).c(modelstruct.order(seq).o(find(resSeq==diff(i),1,'first')))=num2cell(mean(modelstruct.coords(seq).c(resSeq==diff(i),:),1),2);
        end
    end
    modelstruct.resSeq=resSeq;
    modelstruct.centers(seq).c(modelstruct.order(seq).o(ca))=num2cell(modelstruct.coords(seq).c(ca,:),2);
    
    % How many residues lack coordinates
    num_coords=size([modelstruct.centers(seq).c{:}],2)/3;
    if num_coords < modelstruct.numRes(seq).nr
        fprintf('  %s: There are %d residues without coordinates in chain %c.\n',...
            filename, modelstruct.numRes(seq).nr-num_coords, chainID);
    end 
   
    % Save the original numbering
    modelstruct.origSerNo(seq).o=cell((size(modelstruct.centers(seq).c)));
    modelstruct.origSerNo(seq).o(~cellfun(@isempty,modelstruct.centers(seq).c))=num2cell(unique(resSeq));    
    % Save the original residues in letter format
    modelstruct.origRes(seq).r=origSeq;

    % Pairwise distances 
    distances = pdist(modelstruct.coords(seq).c);
    distances2=squareform(distances);
    
    % Find contacting atoms (with the given thershold)
    neighbours = distances<=threshold;
    [r,c] = find(squareform(neighbours));
       
    % Count neighbours of each residue with sparse
    modelstruct.contacts(seq).i =sparse(modelstruct.numRes(seq).nr,modelstruct.numRes(seq).nr);
    modelstruct.contacts(seq).i(1:max(modelstruct.order(seq).o(r)),1:max(modelstruct.order(seq).o(r))) = ...
        sparse(modelstruct.order(seq).o(r),modelstruct.order(seq).o(c),1);
    
    % Determine shortest distances between residues (as triangle matrix)
    modelstruct.residueDistances(seq).d=zeros(modelstruct.numRes(seq).nr,modelstruct.numRes(seq).nr);
    for i=1:modelstruct.numRes(seq).nr
        ord_i=modelstruct.order(seq).o==i; 
        if sum(ord_i)>0
            for j=i+1:modelstruct.numRes(seq).nr
                if sum(modelstruct.order(seq).o==j)>0
                    modelstruct.residueDistances(seq).d(i,j)=min(min(distances2(ord_i,modelstruct.order(seq).o==j)));
                end
            end
        end
    end
    
    % Make the matrix symmetrical
    modelstruct.residueDistances(seq).d=...
        sparse(modelstruct.residueDistances(seq).d+modelstruct.residueDistances(seq).d');
end

%% Find contacts between the the chains
if numSeqs>1
    for i=1:numSeqs-1
        for j=i+1:numSeqs
            distances2=pdist2(modelstruct.coords(i).c, modelstruct.coords(j).c);

            % fnd neighbours with given threshold
            neighbours2 = distances2<=threshold;
            [r,c] = find(neighbours2);

            if isempty(r)
                modelstruct.contacts2(i,j).i=[];
            else
                % use sparse to count how many of each residue are close to others
                modelstruct.contacts2(i,j).i = sparse(modelstruct.numRes(i).nr,modelstruct.numRes(j).nr);
                modelstruct.contacts2(i,j).i(1:max(modelstruct.order(i).o(r)),1:max(modelstruct.order(j).o(c))) = ...
                    sparse(modelstruct.order(i).o(r),modelstruct.order(j).o(c),1);
            end
            
           % Distances between residues
            modelstruct.residueDistances2(i,j).d=...
                zeros(modelstruct.numRes(i).nr,modelstruct.numRes(j).nr);
            for k=1:modelstruct.numRes(i).nr
                ord_k=modelstruct.order(i).o==k; 
                if sum(ord_k)>0
                    for m=1:modelstruct.numRes(j).nr
                        if sum(modelstruct.order(j).o==m)>0
                            modelstruct.residueDistances2(i,j).d(k,m)=...
                                min(min(distances2(ord_k,modelstruct.order(j).o==m)));
                        end
                    end
                end
            end    
        end
    end
end
%% Put all contacts together

if(numSeqs>1)
    
    % starts and ends of each chain after shifting
    starts=struct;
    ends=struct;
    % calculate total number of residues (in all chains)
    tot_res=0;
    for seq=1:numSeqs
        tot_res=tot_res+modelstruct.numRes(seq).nr;
    end
    
    for seq=2:numSeqs
        modelstruct.origSerNo(1).o=[modelstruct.origSerNo(1).o,modelstruct.origSerNo(seq).o];
    end
    
    % contacts within chain in the same matrix, residues are shifted
    % accordingly
    res_contacts=sparse(tot_res,tot_res);
    res_distances=sparse(tot_res,tot_res);
    res_centers=cell(tot_res,1);
    secondary=cell(tot_res,1);
    ind_start=1;
    chains=cell(tot_res,1);
    resOrig_all=[];
    for seq=1:numSeqs
        starts(seq).s=ind_start;
        ind_end=ind_start+modelstruct.numRes(seq).nr-1;
        ends(seq).e=ind_end;
        chains(ind_start:ind_end)=cellstr(char(seq+'A'-1));
        res_contacts(ind_start:ind_end,ind_start:ind_end)= ...
            modelstruct.contacts(seq).i(1:modelstruct.numRes(seq).nr,1:modelstruct.numRes(seq).nr);
        res_distances(ind_start:ind_end,ind_start:ind_end)= ...
            modelstruct.residueDistances(seq).d(1:modelstruct.numRes(seq).nr,1:modelstruct.numRes(seq).nr);
        res_centers(ind_start:ind_end)=modelstruct.centers(seq).c;
        res_all(ind_start:ind_end)=modelstruct.res(seq).r;
        resOrig_all=[resOrig_all,modelstruct.resOrig(seq).r];
        secondary(ind_start:ind_end)=modelstruct.secondary(seq).s;
        ind_start=ind_end+1;
    end

    % contacts between chains
    for i=1:numSeqs-1
        for j=i+1:numSeqs
            if ~isempty(modelstruct.contacts2(i,j).i)
                res_contacts(starts(i).s:ends(i).e,starts(j).s:ends(j).e)= ...
                    modelstruct.contacts2(i,j).i(1:modelstruct.numRes(i).nr,1:modelstruct.numRes(j).nr);
                res_contacts(starts(j).s:ends(j).e,starts(i).s:ends(i).e)= ...
                    modelstruct.contacts2(i,j).i(1:modelstruct.numRes(i).nr,1:modelstruct.numRes(j).nr)';
            end
                res_distances(starts(i).s:ends(i).e,starts(j).s:ends(j).e)= ...
                    modelstruct.residueDistances2(i,j).d(1:modelstruct.numRes(i).nr,1:modelstruct.numRes(j).nr);
                res_distances(starts(j).s:ends(j).e,starts(i).s:ends(i).e)= ...
                    modelstruct.residueDistances2(i,j).d(1:modelstruct.numRes(i).nr,1:modelstruct.numRes(j).nr)';
        end
    end 
else
    tot_res=modelstruct.numRes(seq).nr;  
    res_contacts=modelstruct.contacts(seq).i(1:tot_res,1:tot_res);
    res_distances=modelstruct.residueDistances(seq).d;
    res_centers=modelstruct.centers(seq).c;
    res_all=modelstruct.res(seq).r;
    resOrig_all=modelstruct.origRes(seq).r;
    starts.s=1;
    ends.e=tot_res;
    secondary=modelstruct.secondary(seq).s;
    chains=cell(tot_res,1);
    [chains{:}]=deal(chainID);
    modelstruct.origSerNo(1).o=modelstruct.origSerNo(seq).o;
end

%% everything in good order
structure.numResidues=tot_res;
res_contacts(logical(eye(size(res_contacts))))=0; % residues should not be in contact with themselves
structure.contacts=full(res_contacts);
structure.distances=full(res_distances);
structure.sequence(tot_res)=struct('residues',[],'origSerNo',[],'letters',[],'chain',[],'secondary',[],'centers',[]);
[structure.sequence(:).residues]=res_all{:};
structure.origSequence=resOrig_all;
temp=AA2letters(res_all,true);
[structure.sequence(:).letters]=temp{:};
[structure.sequence(:).chain]=chains{:};
[structure.sequence(:).secondary]=secondary{:};
[structure.sequence(:).centers]=res_centers{:};
[structure.sequence(:).origSerNo]=modelstruct.origSerNo(1).o{:};
structure.chains(numSeqs)=struct('starts',[],'ends',[]);
[structure.chains(:).starts]=starts(:).s;        
[structure.chains(:).ends]=ends(:).e;
structure.score=[modelstruct.score(:).s];

% Add amino acid properties
structure=AAProperties(structure);

fprintf('\n')
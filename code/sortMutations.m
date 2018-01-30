function [mutationsSorted, ddGSorted, indexes, missing]=sortMutations(mutations,ddGs, ...
    WTSeq,removeWTs,removeSelfMuts,removeNoMatch,checkMissing,allowMissingSelfMuts,verbal)
%sortMutations by position number and then by mutated residue. 
% Organise multiple simultanious mutations in increasing order by the
% position numbers. If mutations defined with PDB numbering (that does not
% correspond to the sequence numbering) is used, the WTSeq should not be
% given and cannot check if mutations match the WTseq.
%
%   INPUT:  mutations = cellstr of mutations.
%           ddGs = ddGs corresponding to mutations.
%           WTSeq = wild type sequence. Optional, needed if one wants to
%              check that mutations are correct and to check for missing
%              mutations.
%           removeWTs = true if wild type mutations should be removed.
%           removeSelfmuts = true if mutations such as T53T or A1A should
%               be removed.
%           removeNoMatch = true, if mutations that do not match the WT
%               sequence should be removed. Otherwise warning is given.
%           checkMissing = check if all possible mutations are done, list
%               missing mutations in missing.
%           allowMissingSelfMuts = true if when checking for missing
%               mutations self mutations should not be considered as
%               missing.
%   
%   OUTPUT: mutationsSorted = sorted mutations
%           ddGSorted = ddGs in the same order as the sorted mutations
%           indexes = indexes of the orginal mutations that corrrespond to
%               the sorted mutations.
%           missing = list of missing mutations.

if nargin<7 || isempty(checkMissing)
    checkMissing=false;
end
if nargin<9
    verbal=false;
end
if removeWTs
    remove=cellfun(@(x) contains(x,'wild'),mutations) | cellfun(@(x) contains(x,'WILD'),mutations);
    if sum(remove)>0
        disp('The following mutations were considered as WTs and removed:')
        ind=find(remove,1,'last');
        tmp=remove;tmp(ind)=false;
        tmp2=cellfun(@(x) [x,', '],mutations(tmp),'un',0);
        disp([tmp2{:},mutations{ind}])

    end
else
    remove=false(length(ddGs),1);
end

n=numel(mutations);
sortThis=zeros(n,2);
if isempty(WTSeq)
    if checkMissing
        error('WT sequence is needed to check missing mutations');
    end
    isSeq=false;
else
    isSeq=true;
    if isletter(WTSeq)
        WTSeq=aa2int(WTSeq);
    end
end

inds_m=find(~remove);
for i_m=1:length(inds_m)
    m=inds_m(i_m);
    mut=mutations{m};
    mut=mut(mut~=' '& mut~=',');
    % check that the mutations are valid. Take the residue numbers
    letter_locs=isletter(mut);
    letter_inds=find(letter_locs);

    % handle wild type
    if contains(mut,'WILD')|| contains(mut,'wild')
        warning('mutation %s at index %d is considered as the WT sequence',mut,m)
    % handle mutations with special characters
    elseif any(~ismember(mut,'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'))
        indpdbs=strfind(mut,'(PDB:');
        if ~isempty(indpdbs)
            indpdbs=indpdbs+5;
            newmut='';
            for i=1:length(indpdbs)
                iis=indpdbs(i)-1+find(letter_locs(indpdbs(i):end),2,'first');
                newmut=[newmut, mut((iis(1):iis(2)))];
            end
            warning('Changed %s to %s',mut,newmut);
            mut=newmut;
            letter_locs=isletter(mut);
            letter_inds=find(letter_locs);
        else
            ind=find(~ismember(mut,'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'),1,'first');
            mut_old=mut;
            mut=mut(1:ind-1);
            letter_locs=isletter(mut);
            letter_inds=find(letter_locs);

            if mod(sum(letter_locs),2)==0
                warning('Mutation %s is invalid. Changed to %s',mut_old,mut)
            else
                error('Invalid mutation: %s.',mut)
            end
        end
    % handle mutations with extra letters
    elseif mod(sum(letter_locs),2)~=0 && length(letter_locs)>sum(letter_locs)
        mutparts=strsplit(mutations{m},{' ',','});
        if ~all(cellfun(@(x) sum(isletter(x))==0||sum(isletter(x))==length(x),mutparts))
            while mod(sum(letter_locs),2)~=0
                % if two last are letters
                if letter_locs(end-1)
                     inds=find(letter_locs,3,'last');
                     if double(mut(inds(2)))<double('Q') % is Q good? 
                         newNum=num2str(str2double(mut(inds(1)+1:inds(2)-1))...
                         +(double(mut(inds(2)))-double('A')+1)/100);
                     else
                         newNum=num2str(str2double(mut(inds(1)+1:inds(2)-1))...
                             +(double(mut(inds(2)))-double('A')+1)/100-0.999);
                     end
                     mut=[mut(1:inds(1)),newNum,mut(inds(3))];
                else 
                    locX=strfind(letter_locs, [1 1 1]);
                    for i=1:length(locX)
                        nloc=find(letter_locs(1:locX(i)-1),1,'last')+1;
                        if double(mut(locX(i)))<double('Q') % is Q good? 
                             newNum=num2str(str2double(mut(nloc:locX(i)-1))...
                             +(double(mut(locX(i)))-double('A')+1)/100);
                        else
                             newNum=num2str(str2double(mut(nloc:locX(i)-1))...
                                 +(double(mut(locX(i)))-double('A')+1)/100-0.999);
                        end
                        mut=[mut(1:nloc-1),newNum,mut(locX(i)+1:end)];
                    end
                end
                letter_locs=isletter(mut);
                letter_inds=find(letter_locs);
            end
        else
            warning('Additions are not allowed. Removed mutation %s',mutations{m});
            remove(m)=true;
        end
    end
    
    if remove(m)
        sortThis(m,:)=Inf;
    elseif ~any(ismember(mut,'0123456789'))
        sortThis(m,1)=10e6;
        sortThis(m,2)=10e6;
        mutations{m}=mut;
    else
        n_muts=sum(letter_locs)/2;
        if mod(sum(letter_locs),2)~=0
            disp(mut)
        end
        
        % Check for self mutations
        if removeSelfMuts
            l_is=find(letter_locs);
            for i=1:2:2*n_muts-1
                if mut(l_is(i))==mut(l_is(i+1))
                    remove(m)=true;
                    sortThis(m,:)=Inf;
                    if verbal; fprintf('Removed mutation %s (self-mutation not allowed)\n',mut);end
                    break;
                end
            end
        end
        
        if ~remove(m)
            mut_m=cell(n_muts,1);
            nums=zeros(n_muts,1);
            for i=1:2:n_muts*2
                num=str2double(mut(letter_inds(i)+1:letter_inds(i+1)-1));
                if isSeq && WTSeq(num)~=aa2int(mut(letter_inds(i)))
                    if removeNoMatch
                        remove(m)=true;
                        sortThis(m,:)=Inf;
                        if verbal;fprintf('Removed mutation %s (does not match the WT, the original AA is %c)\n',mut, int2aa(WTSeq(num)));end;
                        break;
                    else
                        warning('Mutation %s does not match the WT, the original AA is %c',mut, int2aa(WTSeq(num)));
                    end
                end
                mut_m{floor(i/2)+1}=[mut(letter_inds(i)) num2str(num) mut(letter_inds(i+1))];
                nums(floor(i/2)+1)=num;
            end

            if ~remove(m)
                [nums,ii]=sort(nums);
                mut_m=mut_m(ii);
                mut=[mut_m{:}];
                sortThis(m,1)=nums(1);
                sortThis(m,2)=aa2int(mut_m{1}(end));
                mutations{m}=mut;
            end
        end
    end
end

sortThis(remove,:)=Inf;
[sortThis,I]=sortrows(sortThis,[1,2]);

indexes=I(1:end-sum(remove));
mutationsSorted=mutations(indexes);
ddGSorted=ddGs(indexes);

if checkMissing
    missing={};
    count=0;
    prevA=20;
    prevR=0;
    m=1;
    if verbal;disp('Missing mutations:');end
    while m<=n

        % If we would expect to see a self mutation, adjust prevA and prevR accordingly
        if allowMissingSelfMuts
            if prevA==19 && WTSeq(prevR)==20 && WTSeq(prevR+1)==1 % both V?V and A?A
                prevA=1;prevR=prevR+1;
            elseif prevA<20 && WTSeq(prevR)==prevA+1 % any one but A?A
                prevA=prevA+1;
            elseif prevA==20 && WTSeq(prevR+1)==1 % just A?A
                prevA=1; prevR=prevR+1;
            end 
        end
        % Basic case next mutation or two of the same mutations
        if (sortThis(m,1)==prevR && (sortThis(m,2)==(prevA+1)|| sortThis(m,2)==prevA))... % same residues
                || (sortThis(m,1)==(prevR+1) && sortThis(m,2)==1 && prevA==20) % consecutive residues
            prevR=sortThis(m,1);
            prevA=sortThis(m,2);
            m=m+1;
        % Something is missing
        else
            if prevA==20
                prevA=1;
                prevR=prevR+1;
            else
                prevA=prevA+1;
            end
            count=count+1;
            missing{count,1}=[int2aa(WTSeq(prevR)) num2str(prevR) int2aa(prevA)];
            if verbal; disp(['  mutation: ' missing{count,1}]); end
        end
    end
    if ~(sortThis(m-1,1)==length(WTSeq) && (sortThis(m-1,2)==20 || ...
            (WTSeq(sortThis(m-1,1))==20 && sortThis(m-1,2)==19 && allowMissingSelfMuts)))
        if verbal;disp(['  Everything after mutation ' mutationsSorted{m-1} ' is missing']);end
        % add to missing
        for mm=sortThis(m-1)+1:length(WTSeq)
            AA_mm=WTSeq(mm);
            rangeAA=1:20;
            if allowMissingSelfMuts; rangeAA(AA_mm)=[];end
            beginning=[int2aa(AA_mm),num2str(mm)];
            for j=rangeAA
                count=count+1;
                missing{count,1}=[beginning, int2aa(j)];
            end
        end
    end
    
    if isempty(missing)
        if verbal;disp('  none');end
    end
end


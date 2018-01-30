function graph=createGraph(protein,label_types)
%% CREATEGRAPH creates a graph-structure of given proteins
% that is in the correct form for graphkernels
%
% INPUT: structure = structure containing the proteins
%        label_types = labels to give for each node (charge, size, hydro,
%        cons), optional
%
% OUTPUT: graphs = structure array with the graph on every row
% 
graph=struct;

vecVals=false;
if nargin>1
    vecVals=true;
    num_types=numel(label_types);
    type_locs=zeros(num_types,1); % column numbers in protein.all
    for i=1:num_types
        if strcmp(label_types{i},'charge')
            type_locs(i)=6; 
        elseif strcmp(label_types{i},'size')
            type_locs(i)=5;
        elseif strcmp(label_types{i},'hydro')
            type_locs(i)=7;
        elseif strcmp(label_types{i},'cons')
            type_locs(i)=4;
        end
    end
end

N=size(protein.contacts,1);

% Adjacency matrix
graph.am=double(protein.contacts>0);

% Node labels: amino acids as integers
graph.nl.values=double(aa2int([protein.sequence(1:N).letters]'));

% Edges and their labels (just ones atm)
[j,k,v]=find(sparse(graph.am));
graph.el.values=[j,k,v];

% Adjacency list: neighbours of every node
graph.al=cell(N,1);
for l=1:N
    graph.al{l}=find(graph.am(l,:));
end

% Add requested vecvalues
if vecVals
    graph.nl.vecvalues=protein.all(:,type_locs);
    graph.nl.vecvalues(:,type_locs==5)=(graph(i).nl.vecvalues(:,type_locs==5)-48)/115;
    graph.nl.vecvalues(:,type_locs==6)=(graph(i).nl.vecvalues(:,type_locs==6)+1)/2;
    graph.nl.vecvalues(:,type_locs==7)=(graph(i).nl.vecvalues(:,type_locs==7)+1)/2;
end
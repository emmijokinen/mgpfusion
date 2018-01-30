%% This is an example for using mGPfusion

% set path for code
addpath(genpath('code'));

% load stability data
load data/stability/ddg_protherm.mat 
load data/stability/ddg_rosetta_single.mat 
load data/stability/ddg_rosetta_multi.mat

% selected substitution matrices
load data/subMats/subMats.mat
accessNos=subMats(:,2);
descsriptions=subMats(:,3);
subMats=subMats(:,1);
% subMats={normBLOSUM62()}; % normalised blosum62 matrix for 20 naturally occurring amino acids
M=length(subMats);

% create data structure, where all necessary data for selected proteins
% will be gathered
pdbs=ddg_protherm(:,1);
data=struct('pdb',pdbs);
% select chain for each protein (to be used from the PDB structure)
% e.g. 1 for first, 2 for second chain and '1', '2', 'A', for a specific name
chains={1,'A',1,1,'I','A',1,1,1,1,1,'A','A','A',1};

%% add structure data
% handle PDB of a protein and gather information
for i=1:length(data)
    pdb=pdbs{i};
    % If you have a pdb-file placed in data/pdb, you can also use e.g.
    %pdb=['data/pdb/' data(i).pdb '.pdb'];
    cutoff=5; % threshold distance for contact (Ångström)
    chain=chains{i}; % amino acid chain of the PDB to consider (letter or number, 0 for all chains)
    prot=protein(pdb,cutoff,chain,{[]});

    % select relevant information from prot
    data(i).prot=cell(1,4);
    data(i).prot{1,1}=prot;
    data(i).prot{1,2}=[prot.sequence(:).letters];

    % get graph form of the protein
    graph=createGraph(prot);
    data(i).prot{1,3}=graph.al; % adjacency list, i.e. neighbours for each residue
    data(i).prot{1,4}=chain;
end
%% add stability data
for i=1:length(data)
    n=length(ddg_protherm{i,2});
    data(i).numMuts=n;
    data(i).mutations=cell(n,3);
    
    mutE=ddg_protherm{i,2}(:,1);
    mutsR=ddg_rosetta_single{i,2}(:,1);
       
    data(i).mutations(:,1)=mutE; 
    data(i).mutations(:,2)=ddg_protherm{i,2}(:,2); % experimental ddg-values
    % corresponding rosetta ddg-values with default scaling
    I=matchSingleMuts(mutE,mutsR);
    data(i).mutations(I>0,3)=ddg_rosetta_single{i,2}(I(I>0),2);
    
    data(i).rosetta_single=ddg_rosetta_single{i,2};
    data(i).rosetta_multi=ddg_rosetta_multi{i,2};
end
 
%% Train model
f=@(a,b,c,d,x) b*x + a.*exp(c*x) + d; % function for transformation

i=6; % index of protein in data
    mutE=data(i).mutations(:,1); % mutations with experimentally measured stability values
    yE=[data(i).mutations{:,2}]'; 
    NE=length(yE)+1; % number of experimental mutations + wild type
    mutR=data(i).rosetta_single(:,1); % mutations with simulated stability values
    yR=[data(i).rosetta_single{:,2}]';
    NR=length(yR); % number of simulated mutations
    N=NE+NR; % number of all mutations
    
    seq=data(i).prot{1,2}; % amino acid sequence of the protein
    pdb=data(i).pdb; % PDB id of the protein

    % Calculate kernel matrix
    al=data(i).prot{1,3}; % adjacency list of protein's residues
    %K=mWDK(seq,[mutE;mutR],al,1,1,B0); % WT is added in the beginning 
    Ks = cellfun(@(S) mWDK(seq,[mutE;mutR],al,1,1,S),subMats,'un',0);
    
    % select experimental data for training and testing
    trE=true(NE-1,1);
    trE(randperm(NE-1,20))=false; % select 20 data points for testing 
    tsE=~trE;
    % use WT, selected experimental, and all simulated data for training
    tr=[true; trE; true(NR,1)]; ts=~tr;
    
    % Select all mutations from training data with experimental and
    % simulated stability values
    I=trE & ~cellfun(@isempty, data(i).mutations(:,3));
    yEm=yE(I);
    yRm=[data(i).mutations{I,3}]';

    % Use Metropolis-Hastings to sample parameters for the transformations
    % of the simulated data. Also get standard deviation of the parameters
    [pars,stdi]=sample_params(yRm,yEm,yR);    

    yS=f(pars(1),pars(2),pars(3),pars(4),yR);
    y=[0;yE;yS];
    
    % normalise y
    my=mean(y(tr));
    yn=y-my;
    ymax=max(abs(yn(tr)));
    yn=yn/ymax;

    Ktr=cellfun(@(K) K(tr,tr),Ks,'un',0);
    
    nWT=10e-6; % standard deviation for the error of wild type ddg
       
    x0=[0.075/ymax;ones(M,1)/M*0.9;zeros(M,1);ones(M,1);0.1/ymax;1.1];
    model=gpmkl([0;yE(trE)],yS,Ktr,'nwrt','tmp',x0,stdi/ymax,my,ymax,nWT);
    model.theta=pars; % parameters for transformation
    model.stdT=stdi; % standard deviation for transformation
    model.ymax=ymax; % parameters for normalising y
    model.my=my;
    model.x0=x0; % used starting values for optimised parameters
    
%% Predict

noise=[nWT; model.n*ones(NE-1,1); model.t*model.stdT/model.ymax+model.n+model.r];

Ky=diag(noise.^2);
for m=1:M
    Ky=Ky+model.w(m)*(Ks{m}+model.c(m)).^model.g(m);
end

[mu,stdPred,~,~,~,~]=gppredict(yn,Ky,noise.^2,tr,ts,0);
ypred=mu*ymax+my;
model.stdPred=stdPred*ymax;

%% Analyze

figure();
scatter(yE(tsE),ypred)
rho=corr(yE(tsE),ypred); % correlation
rmse=sqrt(immse(yE(tsE),ypred)); % root-mean-square error
title([pdb ': $\rho$=' num2str(rho,4) ', rmse=' num2str(rmse,4)],'Interpreter','latex')
set(gca,'DataAspectRatio',[1 1 1])
lsline


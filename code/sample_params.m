function [pars, stdi, smpls] = sample_params(yRm, yEm, yR, NS, burnin)
%% SAMPLE_PARAMS uses Metropolis-Hastings to sample parameters a,b,c,d 
% for scaling of Rosetta data with g(x)= b*x + a*e^(c*x) + d
%
% INPUT:    yRm = Rosetta ddGs that match to yEm
%           yEm = experimentally measured ddGs that match to yRm
%           yR  = all Rosetta ddGs (can be empty if stdi is not needed)
%           
% OUTPUT:   pars = parameters [a b c d]
%           stdi = standard deviation of transformations at each point of yR
%           smpls = all sampled parameters (NS x 4)

if nargin<4 || isempty(NS); NS=10000; end
if nargin<5 || isempty(burnin); burnin=500; end
start=[0.5 1.5 0.2 -0.9];

% Function
f=@(a,b,c,d,x) b*x+a.*exp(c*x)+d; %either multiple params or multiple x possible

% Priors
priorA= @(a) gampdf(a,2,1.5); 
priorB= @(b) betapdf(b/2,1.3,2);
priorC= @(c) betapdf(3.33*c,2,5);
priorD= @(d,a) normpdf(d,-a,0.15);

sigma=0.5;
mypdf=@(par) prod(normpdf(f(par(1),par(2),par(3),par(4),yRm),yEm,sigma))...
    * priorA(par(1)) * priorB(par(2)) * priorC(par(3)) *priorD(par(4),par(1));
proprnd=@(x) x+repmat([.4 .04 .04 .4],size(x,1),1).* rand(size(x))-repmat([.2 .02 .02 .2],size(x,1),1);
try
    % try with default startpoint
    smpl=mhsample(start,NS,'pdf',mypdf,'proprnd',proprnd,'burnin',burnin,'symmetric',1);     
catch
    % try different startpoints if necessary
    fprintf('Starting point was not feasible. Trying different starting points.\n');
    starts=combvec(0:0.5:8,0:0.25:3,0:0.05:0.3,0:-0.5:-8)';
    found=false;
    for i=1:size(starts,1)
        yy=mypdf(starts(i,:));
        if yy>0
            fprintf('New starting point found.\n')
            smpl=mhsample(starts(i,:),NS,'pdf',mypdf,'proprnd',proprnd,'burnin',burnin,'symmetric',1);
            found=true;
            break;
        end
    end
%     i=0;
%     while ~found && i<10^6
%         start=[rand(1)*10, rand(1)*3.5, rand(1)*3.33, rand(1)* (-10)];
%         yy=mypdf(start);
%         if yy>0
%             fprintf('New starting point found.\n')
%             smpl=mhsample(start,NS,'pdf',mypdf,'proprnd',proprnd,'burnin',burnin,'symmetric',1);
%             found=true;
%         end
%         i=i+1;
%         if mod(i,10000)==0;disp(i);end
%     end
%     
    if ~found
        error('Could not find feasible starting point');
    end
end
pars=mean(smpl);

% Calculate standard deviation for yR
if nargout>1
    NR=length(yR);
    yR_i=f(pars(1),pars(2),pars(3),pars(4),yR);
    vari=zeros(NR,1);
    for i=1:NR
        tmp=f(smpl(:,1),smpl(:,2),smpl(:,3),smpl(:,4),yR(i));
        vari(i)=sum((tmp-yR_i(i)).^2);
    end    
    stdi=sqrt(vari/NS);
end

if nargout>2
    smpls=smpl;
end
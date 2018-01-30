function mdl = gpmkl(yE,yR,K,pars,alg,x0,noiseT,my,ymax,nWT,opt_display,tol,max_iter)
% GP-MLL over polynomial kernel: K = w (K+c).^g + noise^2 * I
%
% Parameters
%  'yE'   : (NE x 1) vector of experimental measurements
%  'yR'   : (NR x 1) vector of Rosetta simulations
%  'pars' : 'nwcgr' (n: std of noise w: weights, c: constants, 
%           g: exponents, r: std of Rosetta, t: scaling parameter for
%           std of noise for transformation
%  'alg'  : {'ip','as','sqp','trr'}
%  'x0'   : starting values. If a parameter is not optimised (i.e. not
%       in pars), it will be left to the value given in x0
%  'noiseT': additional elementwise std for noise for Rosetta data 
%  'my'    : coefficient for centering the data (optional, default is mean)
%  'ymax'  : coefficient for scaling the data (optional, default is max(abs([yE; yR]))
%  'opt_display : {'off','none','iter','iter-detailed','notify','notify-detailed',...
%            'final','final-detailed'}
% 'tol' :    tolerance for objective function
% 'max_iter': max number of iterations
% Output
% 'mdl' : structure with fields for optimised parameter values, and mll
	
	% choose optimisation function
	function x = optim(func)
		opts1 = optimoptions('fmincon','GradObj','on','Display',opt_display,...
            'SpecifyObjectiveGradient',true,'MaxIterations', max_iter,...
            'OptimalityTolerance',tol,'MaxFunEvals',10000,...
            'StepTolerance',1e-10,...%'SubproblemAlgorithm','cg',
            'HessianApproximation','lbfgs','ScaleProblem','obj-and-constr');

		switch alg
			case 'ip'
				opts1.Algorithm = 'interior-point';
                x = fmincon(func, x0(Ip),[],[],[],[],lb(Ip),ub(Ip),[],opts1);
			case 'as'
				opts1.Algorithm = 'active-set';
				x = fmincon(func, x0(Ip),[],[],[],[],lb(Ip),ub(Ip),[],opts1);
			case 'sqp'
				opts1.Algorithm = 'sqp';
				x = fmincon(func, x0(Ip),[],[],[],[],lb(Ip),ub(Ip),[],opts1);
			case 'trr'
				opts1.Algorithm = 'trust-region-reflective';
                opts1.HessianApproximation = 'finite-difference';
				x = fmincon(func, x0(Ip),[],[],[],[],lb(Ip),ub(Ip),[],opts1);
            case 'tmp'
                opts.verbose=0; opts.maxIter=500;
                x = minConf_TMP(func,x0(Ip),lb(Ip),ub(Ip),opts);    
		end
    end

    y=[yE;yR];
    NE = length(yE);
	NR = length(yR);
    if isempty(noiseT); noiseT=zeros(NR,1); end

    if nargin<8||isempty(my); my=mean(y); end
    y = y - my;
    if nargin<9||isempty(ymax); ymax=max(abs(y)); end
    y = y / ymax;   
    if nargin<10||isempty(nWT); nWT=10e-6; end
    if nargin<11||isempty(opt_display)||~ismember(opt_display,{'off','none','iter',...
            'iter-detailed','notify','notify-detailed','final','final-detailed'})
        opt_display='off';
    end
    if nargin<12||isempty(tol); tol=1e-6; end
    if nargin<13||isempty(max_iter); max_iter=1000; end
	
    M=length(K);
    
	% parameters
	if isempty(pars); pars = 'wg'; end
    if isempty(alg);  alg = 'ip'; end
    
    if isempty(x0)                   
        % [n,w,c,g,r,t]
        x0 = [0.05/ymax; 0.6/M*ones(M,1); zeros(M,1); ones(M,1); 0.3/ymax; 1.01];
    end
    
    if ~ismember('g',pars)
        x0(2*M+2:3*M+1)=max(x0(2*M+2:3*M+1),ones(M,1));
    else
        % if g is on the lower bound 1 or below it, set to 1.001
        x0(2*M+2:3*M+1)=max(x0(2*M+2:3*M+1),1.001*ones(M,1));
    end
    
    lb = [0.001/ymax; zeros(M,1); zeros(M,1); ones(M,1); 0.001/ymax; 1];
    ub = [10/ymax; 20*ones(M,1); ones(M,1); 30*ones(M,1); 10/ymax; 100];
      
    % don't precompute log(K) to save memory (and if c is optimised)
    logK=cellfun(@(Ki) log0(Ki+x0(3)),K,'un',0);
    
    % optimization criteria
    func = @(x) mll(x,K,y,pars,NE,noiseT,x0,ymax,nWT,logK);

    % define which parameters will be optimised
    Ip=false(3*M+3,1);
    Ip(1)=ismember('n',pars)'; 
    Ip(2:M+1)=ismember('w',pars)'; 
    Ip(M+2:2*M+1)=ismember('c',pars)';
    Ip(2*M+2:3*M+1)=ismember('g',pars)';
    Ip(3*M+2)=ismember('r',pars)'; 
    Ip(3*M+3)=ismember('t',pars)';
    
	% optimise parameters
	x = optim(func);
    mdl = [];
	[mdl.n,mdl.w,mdl.c,mdl.g0,mdl.r,mdl.t] = extractvalues(x,pars,x0); 
	mdl.g = round(mdl.g0); 
     
    % If weight below 1e-6 have been set to zero during optimisation, set
    % them to zero in final results too for consistency
%     Iw=mdl.w<1e-6;
%     mdl.w(Iw)=0;
%     x=[x(1);mdl.w;x(M+2:end)];
	
    mdl.mll = mll(x,K,y,pars,NE,noiseT,x0,ymax,nWT);
end

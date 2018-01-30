function [f,gr] = mll(x,K,y,pars,NE,noiseT,x0,ymax,nWT,logK)
% compute MLL and its gradient
 
    % Parameters and derivative for noise priors and 
    a_E=2.5; b_E=0.02;   a_R=50; b_R=0.007;
    dloggam= @(x,a,b) (a*b-b-ymax*x)/(b*x);
    
	N = length(y);
	M = length(K);
    
    if ~exist('pars','var')
		pars = 'nwg';
    end
	
    % extract values
    [nE,w,c,g,nR,t] = extractvalues(x,pars,x0);

    %Ip = ismember('nwcgrt',pars)';
    Ip=false(3*M+3,1);
    Ip(1)=ismember('n',pars); 
    Ip(2:M+1)=ismember('w',pars); 
    Ip(M+2:2*M+1)=ismember('c',pars);
    Ip(2*M+2:3*M+1)=ismember('g',pars);
    Ip(3*M+2)=ismember('r',pars); 
    Ip(3*M+3)=ismember('t',pars);
  
    noise=[nWT; nE*ones(NE-1,1); t*noiseT+nE+nR];
    %Ky=w*(K+c).^(g) + diag(noise.^2);
   
    Ky=diag(noise.^2);
    Kgs=cellfun(@(Ki,gi) Ki.^gi,K,num2cell(g),'un',0);
    for m=1:M
        Ky=Ky+w(m)*(Kgs{m});
    end
        
    % compute MLL
    try
        try 
            logdet=2*sum(log(diag(chol(Ky))));
        catch
            [~, U, P] = lu(Ky);
            du = diag(U);
            tmp = det(P) * prod(sign(du));
            logdet = log(tmp) + sum(log(abs(du)));
        end
        f = -( -0.5*y'/Ky*y -0.5*logdet -N/2*log(2*pi) );
        if ismember('n',pars)
            f = f - log(ymax*gampdf(nE*ymax,a_E,b_E));
        end
        if ismember('r',pars)
            f = f - log(ymax*gampdf(nR*ymax,a_R,b_R));
        end
    catch
        f=inf;
    end         
	if nargout == 1; return; end
	
	% compute gradient
	a = Ky\y;
    iKy=inv(Ky);
	A = a*a' - iKy;

    dnE= 0; dnR= 0; dt = 0;
	dw = zeros(M,1);
    dc = zeros(M,1); 
    dg = zeros(M,1); 
    
    diagA=diag(A);
    if ismember('n',pars)  
        dnE = nE*sum(diagA(2:end)) + sum( (nR+t*noiseT) .* diagA(NE+1:end) ) + dloggam(nE,a_E,b_E); 
    end
    if ismember('r', pars)
        dnR = sum((nR +nE+t*noiseT) .* diagA(NE+1:end) ) + dloggam(nR,a_R,b_R);
    end
    if ismember('t',pars)
        dt = sum((t*noiseT.^2+(nE+nR)*noiseT) .* diagA(NE+1:end) );
    end
    if ismember('c',pars)
        %dc = 0.5 * sum(diag(A * w*g*(K+c).^(g-1)));
        dc = 0.5 * cellfun(@(Ki,wi,ci,gi) sum(diag(A *wi*gi*(Ki+ci).^(gi-1))),K,num2cell(w),num2cell(c),num2cell(g));
    end
    if ismember('w',pars) && ismember('g',pars)
        for i=1:M
            dw(i)=0.5*trace(A*Kgs{i});
            d=w(i)*Kgs{i}.*logK{i};
            dg(i)= 0.5* ( a'*d*a - trace(iKy*d) );
            %dg(i)=0.5* ( a'*d*a - trace(Ky\d) ); %sligthly more accurate but slower
        end
    else
        if ismember('w',pars)
            dw = 0.5 * cellfun(@(Kg) trace(A * Kg), Kgs);
        end
        if ismember('g',pars)
            for i=1:M
                d=w(i)*Kgs{i}.*logK{i};
                dg(i)= 0.5* ( a'*d*a - trace(Ky\d) ); % iKy*d -> Ky\d
                %dg(i)=  0.5* yiKy * ( d ) * a - 0.5*trace(iKy*d); % This may be more accurate (order 1e-6)
            end
        end
    end
        
	% minimize
    gr = -[dnE;dw;dc;dg;dnR;dt];
	gr = gr(Ip);
	
end
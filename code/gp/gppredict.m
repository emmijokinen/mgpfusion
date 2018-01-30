function [mu,std,mse,cc,mlpd,mll] = gppredict(y,K, nstd, tr, ts,my)
% GP prediction with full 'K' matrix and full 'y' 

	Ntr = nnz(tr);
	Ktt = K(tr,tr);
	Kst = K(ts,tr);
	Kss = K(ts,ts);
    if numel(nstd)==1
        Ky = Ktt + nstd^2*eye(Ntr);
    else
        Ky = Ktt + diag(nstd(tr));
    end
    
	%zs = zeros(Ntr,1);
    if nargin<6
        my = mean(y(tr));
        disp(my)
    end 
	
	A = Kst / Ky;

	mu = A*(y(tr) - my) + my;
	std = sqrt(diag( Kss - A*Kst' ));
	
	mse = mean((mu - y(ts)).^2);
	cc = corr( mu, y(ts) );
%     if numel(nstd)==1
%         mlpd = logmvnpdf( y(ts), mu, diag(std.^2 + nstd^2) );
%     else
%         mlpd = logmvnpdf( y(ts), mu, diag(std.^2 + nstd(ts).^2) );
%     end
%	mll = logmvnpdf( y(tr)-my, zs, Ky );

mlpd=[]; mll=[];
end

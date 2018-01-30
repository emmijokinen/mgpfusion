function X = log0(X)
	X = log(X);
	X(isinf(X)) = 0;
end
	

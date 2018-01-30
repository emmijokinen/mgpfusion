function [n,w,c,g,r,t] = extractvalues(x, pars,x0)
  
    M=(length(x0)-3)/3;
    Ip=false(3*M+3,1);
    Ip(1)=ismember('n',pars)'; 
    Ip(2:M+1)=ismember('w',pars)'; 
    Ip(M+2:2*M+1)=ismember('c',pars)';
    Ip(2*M+2:3*M+1)=ismember('g',pars)';
    Ip(3*M+2)=ismember('r',pars)'; 
    Ip(3*M+3)=ismember('t',pars)';
    
    x0(Ip)=x;
    [n,w,c,g,r,t]=deal(x0(1),x0(2:M+1),x0(M+2:2*M+1),x0(2*M+2:3*M+1),x0(3*M+2),x0(3*M+3));
end
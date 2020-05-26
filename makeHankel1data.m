function Hdata = makeHankel1data(N_multi,z)


Hdata=zeros(2*N_multi+1,1);

Hdata(N_multi+1+0)=besselh(0,1,z);
Hdata(N_multi+1+1)=besselh(1,1,z);

for n=2:N_multi
   
    Hdata(N_multi+1+n)=-Hdata(N_multi+1+n-2)+2*(n-1)/z*Hdata(N_multi+1+n-1);
    
end

for n=-1:-1:-N_multi
   
    Hdata(N_multi+1+n)=-Hdata(N_multi+1+n+2)+2*(n+1)/z*Hdata(N_multi+1+n+1);
    
end

% -besselj(n-2,x)+2*(n-1)/x*besselj(n-1,x)
% besselj(n,x)
% 2*(n+1)/z*besselj(n+1,z)+besselj(n+2,z)

end
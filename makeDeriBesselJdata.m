function deriJdata = makeDeriBesselJdata(N_multi,z,Jdata)

if N_multi==0

    deriJdata=(besselj(-1,z)-besselj(1,z))/2;
else
    
    
deriJdata=zeros(2*N_multi+1,1);


for n=-(N_multi-1):N_multi-1
   
    deriJdata(N_multi+1+n)=(Jdata(N_multi+1+n-1) - Jdata(N_multi+1+n+1)    )/2;
    
end

n=-N_multi;
deriJdata(1)=   n/z*Jdata(N_multi+1+n)-Jdata(N_multi+1+n+1)  ;

n=N_multi;
deriJdata(2*N_multi+1)= Jdata(N_multi+1+n-1)-n/z*Jdata(N_multi+1+n);

end

end
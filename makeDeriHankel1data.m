function deriHdata = makeDeriHankel1data(N_multi,z,Hdata)

if N_multi==0

    deriHdata=(besselh(-1,1,z)-besselh(1,1,z))/2;
else


deriHdata=zeros(2*N_multi+1,1);


for n=-(N_multi-1):N_multi-1
   
    deriHdata(N_multi+1+n)=(Hdata(N_multi+1+n-1) - Hdata(N_multi+1+n+1)    )/2;
    
end

n=-N_multi;
deriHdata(1)=   n/z*Hdata(N_multi+1+n)-Hdata(N_multi+1+n+1)  ;

n=N_multi;
deriHdata(2*N_multi+1)= Hdata(N_multi+1+n-1)-n/z*Hdata(N_multi+1+n);

end

end
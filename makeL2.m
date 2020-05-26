function L = makeL2(i,j,ci,cj,Ri,Rj,k0,delta,N_multi,SpecialFuncDataLi,SpecialFuncDataLj,Hdata_cicj)



%SpecialFuncDataL=[Jdata_k0R,dJdata_k0R];
Jdata_k0Ri=SpecialFuncDataLi(:,1);
dJdata_k0Ri=SpecialFuncDataLi(:,2);

Jdata_k0Rj=SpecialFuncDataLj(:,1);
dJdata_k0Rj=SpecialFuncDataLj(:,2);


Sk0=zeros(2*N_multi+1);
dSk0=zeros(2*N_multi+1);
zeromatrix=zeros(2*N_multi+1);
const=(-1/2)*1i*pi*Rj;

cicj=cj-ci;
%mag_cicj=norm(cicj);
t_cicj=angle(cicj(1)+1i*cicj(2));

N_multi2 = 2*N_multi+1;

for m=-N_multi:N_multi
for n=-N_multi:N_multi
   
   %Jnk0= besselj(n,k0*R);
   %Jmk0= besselj(m,k0*R);
   %Hnmk0= besselh(n-m,1,k0*mag_cicj);
   %dJmk0= 1/2*(besselj(m-1, k0*R) - besselj(m+1, k0*R));
   Jnk0= Jdata_k0Rj(N_multi+1+n);
   Jmk0= Jdata_k0Ri(N_multi+1+m);
   dJmk0= dJdata_k0Ri(N_multi+1+m);
   
   if j>=i
    Hnmk0= Hdata_cicj(i,j,N_multi2+1+n-m);
   else
    Hnmk0= Hdata_cicj(j,i,N_multi2+1+n-m);
   end
   
   Sk0(m+N_multi+1,n+N_multi+1)=const*Jnk0*Jmk0*(-1)^(n-m)*Hnmk0*exp(1i*(n-m)*t_cicj);
   dSk0(m+N_multi+1,n+N_multi+1)=const*Jnk0*k0*dJmk0*(-1)^(n-m)*Hnmk0*exp(1i*(n-m)*t_cicj);
   
end
end


L=[zeromatrix, -Sk0; zeromatrix, -delta*dSk0];

end
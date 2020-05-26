% make operator M
function M = makeM(R,k0,kb,delta,N_multi,SpecialFuncDataM)


Jdata_k0R=SpecialFuncDataM(:,1);
Jdata_kbR=SpecialFuncDataM(:,2);
Hdata_k0R=SpecialFuncDataM(:,3);
Hdata_kbR=SpecialFuncDataM(:,4);
dJdata_kbR=SpecialFuncDataM(:,5);
dHdata_k0R=SpecialFuncDataM(:,6);



const=(-1/2)*1i*pi*R;

Sk0=zeros(2*N_multi+1);
Skb=zeros(2*N_multi+1);
dSk0=zeros(2*N_multi+1);
dSkb=zeros(2*N_multi+1);

for n=-N_multi:N_multi
   
   Jk0= Jdata_k0R(N_multi+1+n);
   Jkb= Jdata_kbR(N_multi+1+n);
   Hk0= Hdata_k0R(N_multi+1+n);
   Hkb= Hdata_kbR(N_multi+1+n);
   dHk0= dHdata_k0R(N_multi+1+n);
   dJkb= dJdata_kbR(N_multi+1+n);
    
   
   Sk0(n+N_multi+1,n+N_multi+1)=const*Jk0*Hk0;
   dSk0(n+N_multi+1,n+N_multi+1)=const*k0*Jk0*dHk0;
   Skb(n+N_multi+1,n+N_multi+1)=const*Jkb*Hkb;
   dSkb(n+N_multi+1,n+N_multi+1)=const*kb*Hkb*dJkb;
   
end

M=[Skb, -Sk0; dSkb, -delta*dSk0];

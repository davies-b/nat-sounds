function A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy)

N = length(cx);
% R1 = R(1);
% R2 = R(2);

k0=omega*sqrt(rho0/kappa0);
kb=omega*sqrt(rho_b/kappa_b);

Jdata_k0 = [];
Jdata_kb = [];
Hdata_k0 = [];
Hdata_kb = [];
dJdata_k0 = [];
dJdata_kb = [];
dHdata_k0 = [];

for i = 1:N
Jdata_k0 = [Jdata_k0, makeBesselJdata(N_multi, k0*R(i)  )];
Jdata_kb = [Jdata_kb, makeBesselJdata(N_multi, kb*R(i)  )];
Hdata_k0 = [Hdata_k0, makeHankel1data(N_multi, k0*R(i)  )];
Hdata_kb = [Hdata_kb, makeHankel1data(N_multi, kb*R(i)  )];
dJdata_k0 = [dJdata_k0, makeDeriBesselJdata(N_multi,k0*R(i),Jdata_k0(:,i))];
dJdata_kb = [dJdata_kb, makeDeriBesselJdata(N_multi,kb*R(i),Jdata_kb(:,i))];
dHdata_k0 = [dHdata_k0, makeDeriHankel1data(N_multi,k0*R(i),Hdata_k0(:,i))];
end

% dJdata_k0R1=makeDeriBesselJdata(N_multi,k0*R1,Jdata_k0);
%Jdata_k0R=makeBesselJdata(N_multi, k0*R  );
%HdataCiCj

% Jdata_k0R2=makeBesselJdata(N_multi, k0*R2  );
% Jdata_kbR2=makeBesselJdata(N_multi, kb*R2  );
% Hdata_k0R2=makeHankel1data(N_multi, k0*R2  );
% Hdata_kbR2=makeHankel1data(N_multi, kb*R2  );
% dJdata_kbR2=makeDeriBesselJdata(N_multi,kb*R2,Jdata_kbR2);
% dHdata_k0R2=makeDeriHankel1data(N_multi,k0*R2,Hdata_k0R2);
% 
% 
% dJdata_k0R2=makeDeriBesselJdata(N_multi,k0*R2,Jdata_k0R2);

%%
N_bubbles=size(cx,2);
Hdata_cicj= zeros(N_bubbles,N_bubbles,2*(2*N_multi+1)+1);
N_multi2=2*N_multi+1;
for i=1:N_bubbles
   for j=i:N_bubbles
       
       
       cicj=[cx(j)-cx(i),cy(j)-cy(i)];
       mag_cicj=norm(cicj);
       %t_cicj=angle(cicj(1)+1i*cicj(2));

       Hdata_cicj(i,j,N_multi2+1+0) = besselh(0,1,k0*mag_cicj);
       Hdata_cicj(i,j,N_multi2+1+1) = besselh(1,1,k0*mag_cicj);
       
       for n=2:N_multi2
   
            Hdata_cicj(i,j,N_multi2+1+n)=-Hdata_cicj(i,j,N_multi2+1+n-2)+2*(n-1)/(k0*mag_cicj)*Hdata_cicj(i,j,N_multi2+1+n-1);
    
       end
       for n=-1:-1:-N_multi2
   
            Hdata_cicj(i,j,N_multi2+1+n)=-Hdata_cicj(i,j,N_multi2+1+n+2)+2*(n+1)/(k0*mag_cicj)*Hdata_cicj(i,j,N_multi2+1+n+1);
    
       end

       
       %if i ~= j
       %Hdata_cicj(j,i,N_multi2+1+0) = Hdata_cicj(i,j,N_multi2+1+0) ;
       %Hdata_cicj(j,i,N_multi2+1+1) = Hdata_cicj(i,j,N_multi2+1+1) ;
       %end
       
       
       
       
   end
end


%%


SpecialFuncDataM = [];
SpecialFuncDataL = [];
for i = 1:N
SpecialFuncDataM = [SpecialFuncDataM, Jdata_k0(:,i),Jdata_kb(:,i),Hdata_k0(:,i),Hdata_kb(:,i),dJdata_kb(:,i),dHdata_k0(:,i)];
SpecialFuncDataL =[SpecialFuncDataL, Jdata_k0(:,i),dJdata_k0(:,i)];
end


N_oneblock = (2*N_multi+1)*2; 

A=zeros( N_oneblock*N_bubbles );

for i=1:N
    for j=1:N
        if i==j
            A((i-1)*N_oneblock+1:i*N_oneblock,(i-1)*N_oneblock+1:i*N_oneblock)= makeM(R(i),k0,kb,delta,N_multi,SpecialFuncDataM(:,6*(i-1)+1:6*i));
        else
                A((i-1)*N_oneblock+1:i*N_oneblock,(j-1)*N_oneblock+1:j*N_oneblock)= makeL2(i,j,[cx(i),cy(i)],[cx(j),cy(j)],R(j),R(i),k0,delta,N_multi,SpecialFuncDataL(:,2*(i-1)+1:2*i),SpecialFuncDataL(:,2*(j-1)+1:2*j),Hdata_cicj);
        end
    end
end

end
function h = h(t,n,res,V)
V1 = inv(V);
nu = @(n) sum(V1(n,:));
h = zeros(size(t));
I = find(t>0);
h(I) = nu(n)*real(res(n))*exp(imag(res(n))*t(I)).*sin(real(res(n))*t(I));
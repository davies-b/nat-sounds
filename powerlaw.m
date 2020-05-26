function [f,S,p] = powerlaw(a,Fs,k)
% fits a 1/x^gamma relationship to the data a

L = size(a,2);
N = size(a,1);

S_store = zeros(N,(k-1));

for n = 1:N
Y = fft(a(n,:));
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
Sn = abs(P1).^2;

f = Fs*(0:(L/2))/L;

S_store(n,:) = Sn(2:k);
f = f(2:k);
end

S = mean(S_store,1);
p = polyfit(log10(f),log10(S),1);
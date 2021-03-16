close all, clear

%%% Obtain the resonant data: can be generated by calling 
%%% [res,V]=resonances();
load('res22');  
N = length(res);


%%% Obtain the filtered: can be generated by calling, for example
[a,Fs] = audiofilter('samples/trumpet.aif',res,V);
% load('trumpet')  % option to load an example of stored data
L = size(a,2) + 1;
T = 1/Fs;

%% Hilbert transform

a_H = hilbert(a');
a_H = transpose(a_H);

%% Normalisation (to ensure invariance)

lA = log10(abs(a_H));
lA = lA - repmat(mean(lA,2),1,L-1);         % zero mean
A = diag(mean(lA.^2,2).^(-1))*lA;           % unit variance

arg = angle(a_H);
lambda = diff(arg,1,2)/T;
lam = abs(lambda);
lam = lam - repmat(mean(lam,2),1,L-2);      % zero mean
lam = diag(mean(lam.^2,2).^(-1))*lam;       % unit variance

%% Temporal averages of A and lamba
k = Fs/100;     % interval length

avA = zeros(N,L-k-2);
avL = zeros(N,L-k-2);
for i = 1:L-k-2
    avA(:,i) = mean(A(:,i:i+k),2);
    avL(:,i) = mean(lam(:,i:i+k),2);
end
%%
parvec = zeros(6,1);        % vector to store the natural sound parameters
%% Fit PDF to data using least squares
%% p_A

[grid,data,fit,vals] = pdffit_exp(avA(:),1,2);
figure
plot(grid,data)
hold on
plot(grid,fit,'LineWidth',2)
legend('$\langle\log A\rangle$','$p_{A}$','interpreter','latex')
parvec(2:3) = vals(1:2);
set(gca,'TickLabelInterpreter','latex')

%% p_lambda

scale = 1e6;                 % rescale the data for numerical convenience
[grid,data,fit,vals] = pdffit_modcau(scale*avL(:),1,5,1e3);
figure
plot(grid/scale,data*scale)
hold on
plot(grid/scale,fit*scale,'LineWidth',2)
legend('$\langle|\lambda|\rangle$','$p_{\lambda}$','interpreter','latex')
parvec(5:6) = [vals(1)/scale, vals(2)];
set(gca,'TickLabelInterpreter','latex')


%% Fit 1/f power law behaviour
%% A power law
figure
[f,S,p] = powerlaw(A,Fs,2e2);
plot(f,S,'.-');
set(gca,'xscale','log')
set(gca,'yscale','log')
y = polyval(p,log10(f));
hold on
plot(f,10.^y,'LineWidth',2);
parvec(1) = -p(1);
xlabel('$\omega$','interpreter','latex')
ylabel('$S(\omega)$','interpreter','latex')
legend('$\overline{S_{A}}$','$S_{A}$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex')


%% phi power law
figure
[f,S,p] = powerlaw(lam,Fs,1e2);
plot(f,S,'.-');
set(gca,'xscale','log')
set(gca,'yscale','log')
y = polyval(p,log10(f(1,:)));
hold on
plot(f(1,:),10.^y,'LineWidth',2);
parvec(4) = -p(1);
xlabel('$\omega$','interpreter','latex')
ylabel('$S(\omega)$','interpreter','latex')
legend('$\overline{S_{\phi}}$','$S_{\phi}$','interpreter','latex');

function [resonances,V] = resonances()

L = 35e-3;              % length of cochlea
N = 22;                 % number of resonators

s = 1.05;
a = 0.00016;

for i = 1:N
    R(i) = a*s^(i-1);
end

%%% Material parameters
rho0 = 1e3;                 % density of water
kappa0 = 2e9;               % bulk modulus of water
v = sqrt(kappa0/rho0);      % speed of sound in water

rho_b = 1.2;                % density of resonators/air
kappa_b = 1e5;              % bulk modulus of resonators/air
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air

% High contrast parameter \delta
delta=rho_b/rho0;

cx = linspace(0,L,N+2);
cx = cx(2:end-1);
cy = zeros(1,N);

% Maximum order for multipole expansion (n = 0, 1, 2, ..., N_multi)
N_multi = 3;

%% Compute initial guesses for the resonances
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) eigs((MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy)),1,'smallestabs');

x = linspace(1, 2*pi*15000, 800);
init = [];
for correction = [1i 5i 10i 20i]
    y = zeros(1, length(x));
    for i = 1:length(x)
        y(i) = abs(f(x(i)-correction));
    end
    for i = 2:length(x)-1
        if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e0)
            init = [init x(i)-correction];
        end
    end
end

if length(init) < length(R)
    disp('WARNING: fewer than N initial guesses created')
end

init = sort(init);

%% Use Muller's method to find the resonances

distTol = 5e-5; fTol = 1e-5; iterMax = 10;
resonances = [];
n = 1;
for initGuess = init
        
    z0 = initGuess;
    z1 = initGuess - 1i;
    z2 = initGuess - 2i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e0
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi (%.0f Hz) \n'], real(res), imag(res), real(res)/2/pi)
       resonances = [resonances res];
       n = n + 1;
    end
end


%% Computing eigenmodes at centres
%  Note: modes approximtely constant on each resonator, constant is given
%  by corresponding entry of V

if N_multi ~= 3
    disp('Error: set up for N_multi=3 only')
end

N_res = length(resonances);
gridPoints = [cx; cy];                  % Grid for field
V = zeros(N);

for m = 1:N_res
omega = resonances(m);
k = omega/v;                 
kb = omega/v_b;                       

A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);
[U, D] = eig(A);
[~,permutation] = sort(diag(D));
D=D(permutation,permutation); U = U(:,permutation);

phi = [];
psi = [];
N_terms = 2*N_multi+1;

for i = 1:N
    phi = [phi, U((i-1)*14+1:i*14-7,1)];
    psi = [psi, U((i-1)*14+8:i*14,1)];
end

% Calculate field
green = @(k,x,y) -1i/4*besselh(0,1,k*sqrt(x.^2+y.^2));

parfor j = 1:N
    gridPoint = gridPoints(:, j);
    
    % Determine which domain
    I = (gridPoint(1)*ones(1,length(cx))-cx).^2 + (gridPoint(2)*ones(1,length(cy))).^2  <= R.^2 ;
    I = find(I);
    fun1 = @(t) green(kb, gridPoint(1)-cx(I)-R(I).*cos(t), gridPoint(2)-cy(I)-R(I).*sin(t)).*...
        (exp(-3i*t)*phi(1,I)+exp(-2i*t)*phi(2,I)+exp(-1i*t)*phi(3,I)+phi(4,I)+exp(1i*t)*phi(5,I)+exp(2i*t)*phi(6,I)+exp(3i*t)*phi(7,I));
    V(j,m) = integral(fun1, -pi, pi);
end

V = real(V);
end
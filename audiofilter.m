function [a,Fs] = audiofilter(file,res,V)
% returns the sampled signal in the matrix a along with the sampling rate
% Fs, due to a subwavelength scattering transform with resonant properties
% given by res and V

[x1,Fs] = audioread(file);      % get the data and sampling rate for the audio file
x = x1(1:Fs);

T = 1/Fs;                     % Sampling period
L = length(x);                % Length of signal
t = (0:L-1)*T;                % Time vector for signal
N = length(res);              % Number of resonators

a = zeros(N,2*L-1);
for n = 1:N
    % Compute the convolution using the fft/ifft
    a(n,:) = fconv(h(t,n,res,V), x);
end
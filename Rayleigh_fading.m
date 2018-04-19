close all
clear all

% Parameters:
num_of_bits = 100;  % Length of into bits sequence (must be even)
A = 1;              % Amplitude of Tx symbols
sigma_h = 100;      % Variance of Rayleigh channel
sigma_n = 10^-1;    % Variance of AWGN
ISNRdB = 20;        % ISNR at RX per channel

% Simulation:
% TX symbols
TX_constellation_syms = A*[exp(j*pi/4),exp(j*7*pi/4),exp(j*3*pi/4),exp(j*5*pi/4)];
% Trick here! Symbols are ordered to cause Gray coding
TX_bits = round(rand(1,num_of_bits));
tmp1 = reshape(TX_bits, 2, num_of_bits/2);  % Prepare to map bits
tmp2 = tmp1(1,:)*2 + tmp1(2,:) + 1;
TX_syms = TX_constellation_syms(tmp2);      % Gray coding is automatic

% Channel:
h = (randn(1,1) + j*randn(1,1))*sqrt(sigma_h/2);
ISNR = 10.^(ISNRdB/10);
Es = sigma_n * ISNR/(abs(h)^2);

% RX symbols:
n = (randn(1,num_of_bits/2) + j*randn(1,num_of_bits/2)) * sqrt(sigma_n/2);
RX_syms = TX_syms.*h + n;

% Display:
% figure
% polar(angle(TX_syms), abs(TX_syms),'*')
% title('Polar Presentation of TX symbols')
% 
% figure
% subplot(211)
% plot(real(TX_syms),'*')
% grid
% title('Real Part of TX Symbols')
% subplot(212)
% plot(imag(TX_syms),'*')
% grid
% title('Imaginary Part of TX Symbols')

figure
polar(angle(RX_syms), abs(RX_syms),'*')
title('Polar Presentation of RX symbols')

figure
subplot(211)
plot(real(RX_syms),'*')
grid
title('Real Part of RX Symbols')
subplot(212)
plot(imag(RX_syms),'*')
grid
title('Imaginary Part of RX Symbols')

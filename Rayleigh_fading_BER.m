close all
clear all

% Parameters:
SNRdB = [0:2:20];                 % SNR at RX
num_of_bits_per_channel = 10^2;  % Length of into bits sequence (must be even)
channels_per_SNR = 10^3;         % Number of channel instances to simulate per SNR
AWGN_channel = 1;
sigma_h = 1;      % Variance of Rayleigh channel
Es = 1;           % Symbol Energy
SNR = 10.^(SNRdB/10);
BER = zeros(1, length(SNR));

for(i1=1:length(SNR))
    sigma_n = Es*sigma_h / SNR(i1);     % AWGN
    errors = 0;
    for(i2=1:channels_per_SNR)          % This loop can be avoided using Matrix notation
        %% TX symbols
        TX_constellation_syms = sqrt(Es)*[exp(j*pi/4),exp(j*7*pi/4),exp(j*3*pi/4),exp(j*5*pi/4)];
        % Trick here! Symbols are ordered to cause Gray coding
        TX_bits = round(rand(1,num_of_bits_per_channel));
        tmp1 = reshape(TX_bits, 2, num_of_bits_per_channel/2);  % Prepare to map bits
        tmp2 = tmp1(1,:)*2 + tmp1(2,:) + 1;
        TX_syms = TX_constellation_syms(tmp2);      % Gray coding is automatic

        %% Channel:
        if(AWGN_channel == 1)
            h = 1;
        else
            h = (randn(1,1) + j*randn(1,1))*sqrt(sigma_h/2);    % Rayleigh channel
        end 
        % AWGN
        n = (randn(1,num_of_bits_per_channel/2) + j*randn(1,num_of_bits_per_channel/2)) * sqrt(sigma_n/2);
        
        %% RX symbols:
        RX_syms = TX_syms.*h + n;
        
        RX_syms_sync = RX_syms * exp(-j*angle(h));   % Perfect channel estimation
        
        %% ML detection
        b1_est = (sign(real(-RX_syms_sync)) + 1)/2;
        b2_est = (sign(imag(-RX_syms_sync)) + 1)/2;
        RX_bits = reshape([b1_est;b2_est], 1, num_of_bits_per_channel);
        errors = errors + sum(TX_bits~=RX_bits);
    end
    BER(i1) = errors/(num_of_bits_per_channel * channels_per_SNR);
end

    

%Display:
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

figure
semilogy(SNRdB, BER)
hold on

%BER = berfading(SNRdB,'psk',8,1);
BER = berawgn(SNRdB,'psk',8,'nondiff');
semilogy(SNRdB, BER)
legend('Simulation','Theory','Location','Best')
title('Average Bit Error Rate - QPSK in AWGN')
xlabel('SNR[dB]')
ylabel('BER')
grid
axis([0 20 (10^-3) 1])

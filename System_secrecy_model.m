close all
clear all

% Parameters:
num_of_bits = 200;  % Length of into bits sequence (must be even)
A = 1;              % Amplitude of Tx symbols
K = 10;
N = 8;
mu = sqrt( K/(2*(K+1)) );
s = sqrt( 1/(2*(K+1)) );

Ts = 1:100;
yc = 0;
ys = 0;
% Simulation:
% TX symbols
TX_constellation_syms = A*[exp(j*pi/4),exp(j*7*pi/4),exp(j*3*pi/4),exp(j*5*pi/4)];
% Trick here! Symbols are ordered to cause Gray coding
TX_bits = round(rand(1,num_of_bits));
tmp1 = reshape(TX_bits, 2, num_of_bits/2);  % Prepare to map bits
tmp2 = tmp1(1,:)*2 + tmp1(2,:) + 1;
TX_syms = TX_constellation_syms(tmp2);      % Gray coding is automatic

% Receiver Channel:
%Rician fading channel - implementationn using Novel Sum of sinusoids
    theta = -pi + 2.*pi.*rand(N,1);
    theta(1) = pi/4;
    phi = -pi + 2.*pi.*rand(N,1);
    phi(1) = pi/4;
    for n= 1:N
        alpha_n = (2.*pi.*n + theta(n))/N;
        yc = yc + cos(Ts.*cos(alpha_n) + phi(n));
        ys = ys + sin(Ts.*cos(alpha_n) + phi(n));
    end
    Yc = (1./sqrt(N)).*yc;
    Ys = (1./sqrt(N)).*ys;
        
    Zc = (Yc + sqrt(K).*cos(Ts.*cos(theta(1))+phi(1)))./sqrt(1+K);
    Zs = (Ys + sqrt(K).*sin(Ts.*cos(theta(1))+phi(1)))./sqrt(1+K); 
    Z = Zc + j.*Zs;
    
%Compensated symbols
TX_syms_comp = TX_syms .*Z;
% RX symbols:
RX_syms = TX_syms_comp.*exp(-j*angle(Z));


% Channel at eavesdropper
    theta = -pi + 2.*pi.*rand(N,1);
    theta(1) = pi/4;
    phi = -pi + 2.*pi.*rand(N,1);
    phi(1) = pi/4;
    for n= 1:N
        alpha_n = (2.*pi.*n + theta(n))/N;
        yc = yc + cos(Ts.*cos(alpha_n) + phi(n));
        ys = ys + sin(Ts.*cos(alpha_n) + phi(n));
    end
    Yc = (1./sqrt(N)).*yc;
    Ys = (1./sqrt(N)).*ys;
        
    Zc = (Yc + sqrt(K).*cos(Ts.*cos(theta(1))+phi(1)))./sqrt(1+K);
    Zs = (Ys + sqrt(K).*sin(Ts.*cos(theta(1))+phi(1)))./sqrt(1+K); 
    Z_ea = Zc + j.*Zs;
Eavesdropper_syms = TX_syms.*Z_ea;      
Eavesdropper_syms_sync = Eavesdropper_syms .* exp(-j*angle(Z_ea));

%Display:
figure
polar(angle(TX_syms), abs(TX_syms),'*')
title('Polar Presentation of TX symbols')

figure
polar(angle(RX_syms), abs(RX_syms),'*')
title('Polar Presentation of RX symbols')

figure
polar(angle(Eavesdropper_syms), abs(Eavesdropper_syms),'*')
title('Polar Presentation of sybols at Eavesdropper')

figure 
plot(real(RX_syms))
threshold = -0.5 + zeros(100);
hold on;
plot(threshold);
title('Deep fades due to very low Doppler spread')
    xlabel('Time : wd*tau')
    ylabel('Re[RX_symbols]')
    legend('RX signal envelope','threshold','Location','Best')
close all
clear all

%parameters:
num_of_channels = 10^5;
sigma_h = 1;    %variance of Rayleigh channel

%stats
h = (randn(1,num_of_channels)+j*randn(1,num_of_channels))*sqrt(sigma_h/2);
alpha = abs(h);
phi = atan2(imag(h), real(h));

%display:
[y,x] = hist(alpha);
y = y/num_of_channels;
figure
subplot(211)
plot(x,y)
axis([0,4*sigma_h,0,1])
title('Estimated PDF of Channel Magnitude (alpha)')
[y,x] = hist(phi);
y = y/num_of_channels;
subplot(212)
plot(x,y)
axis([-pi,pi,0,1])
title('Estimated PDF of Channel Phase (phi)')
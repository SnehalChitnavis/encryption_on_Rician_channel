close all;
clear all;

% 10MHz
N = 8;
N_theory = inf;

%Initialize
tau = 1000;
%wd = 10;
Tp = 10000;
X = 0:100;
J0=besselj(0,X);
yc = 0;
ys = 0;

K =3;

theta = -pi + 2.*pi.*rand(N,1);
theta(1) = pi/4;

%R_ZZ_th = ((J0 + K.* cos(wd*tau.*cos(theta(1))) + j.*K.*sin(wd*tau.*cos(theta(1))))/(1+K) == 0);
%plot(X,real(R_ZZ_th));
%plot(X,imag(R_ZZ_th));
syms wd;
R_ZZ_th = (J0 + K.* cos(X.*cos(theta(1))))/(1+K);
temp = fzero(@(wd)((besselj(0,wd*tau) + K.* cos(wd*tau.*cos(theta(1))))/(1+K)), 3)

alpha =  temp/(2*pi)
Bd = [0:1000];
Bc = 2*Bd/alpha;
eff = 1 - (2./(alpha.*(Bc./Bd)));
plot((Bd),eff);

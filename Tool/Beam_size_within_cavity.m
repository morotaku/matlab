clear all
close all
clc

%%%Beam parameter(μm)
n=1;
lam=640*10^(-3);

%%%Cavity parameter (μm)
l=8*10^3;
R1=-10^10;
R2=10^10;
f_lens=25.4*10^3

%%Rayleigh length of oscilating mode
z0=sqrt(l*(-R1-l)*(R2-l)*(R2-R1-l)/(R2-R1-2*l)^2)

%%%Beam waist within the cavity
w0=sqrt(lam*z0/(pi*n))

%%% Beam size at position z
z=3;
w=w0*sqrt(1+(z/z0)^2)
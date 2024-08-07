clear all; close all; clc;
R=150;
L=linspace(5,20,100);
A=(-58.*L+29*R+90)./(29.*R);
D=(-58.*L+90)./(29.*R);

plot(L,(A+D)/2)
ylim([-2 2])
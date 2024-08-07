clear all
close all
clc
%initial phase
ini_x=pi;
ini_y=pi;
%phase delay
phi_x=0;
phi_y=-pi/4;
%Intensity ratio
thita=pi/2;
int_x=cos(thita/2);
int_y=sin(thita/2);
for t=0:100
    x=cos(pi/20*t+phi_x+ini_x)*int_x;
    y=cos(pi/20*t+phi_y+ini_y)*int_y;
    quiver(0,0,x,0)
    xlim([-1 1])
    ylim([-1 1])
    hold on
    quiver(0,0,0,y)
    hold on
    quiver(0,0,x,y)
    hold off
    pause(0.1)
end
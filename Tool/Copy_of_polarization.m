clear all
close all
clc
x=linspace(-2,2,100);
y=linspace(-2,2,100);
z=linspace(0,6*pi,100);
Z=linspace(-1,6*pi,100);
ex=1;ey=1;d=pi/2;
ax=ex.*cos(z);
ay=ey.*cos(z+d);
az=sqrt(ax.^2+ay.^2);


plot3(ax,z,ay.*0,'LineWidth',1,'LineStyle','--')
hold on
plot3(ax.*0,z,ay,'LineWidth',1,'LineStyle','--')
hold on
plot3(ax,z,ay,'LineWidth',2,'Color','Green')
hold on
plot3(x./1.5,z.*0+6*pi,ay.*0,'LineWidth',1,'Color','black')
hold on
plot3(ax.*0,z.*0+6*pi,y./1.5,'LineWidth',1,'Color','black')
hold on
plot3(x.*0,Z.*1.25,y.*0,'LineWidth',1,'Color','black')
hold on
xlim([-2 2]);zlim([-2 2]);
axis off

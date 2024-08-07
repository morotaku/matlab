clear all
close all
clc
%Radial index
l1=5;
l2=-l1;
l3=0;
l4=-l3;
l5=0;
l6=-l5;
l7=35;
l8=-l7;
l9=39;
l10=-l9;
%l=[l1,l2,l3,l4,l5,l6,l7,l8,l9,l10];

%Azimuth index
p1=1;
p2=0;
p3=2;
p4=3;
p5=4;
p6=3;
p7=3;
p8=3;
p9=3;
p10=3;
%p=[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];

%Parameters
w=1; %Beam waist
lam=0.64; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
z=input('position: '); %Beam position 
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
%Z=z/zr

%x-y　coordinate
N=1000;
L=10; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
%Polar coodinate
%r=sqrt(x^2+y^2);
%phi=atan2(y,x);
[phi,r] = cart2pol(x,y);

p=0;
l1=-2;
l2=1;
C=(2*factorial(p)/(pi*factorial(p+abs(l1))))^0.5;
Gouy=-1j*((-(l1*phi)-((2*p+abs(l1)+1)*atan(z/zr)))+k*z);%+k*(r.^2/(2*R))+k*z);
LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l1).*laguerre(p, l1, 2*r.^2/W^2).*exp(Gouy);


C=(2*factorial(p)/(pi*factorial(p+abs(l2))))^0.5;
Gouy=-1j*((-(l2*phi)-((2*p+abs(l2)+1)*atan(z/zr)))+k*z);%+k*(r.^2/(2*R))+k*z);
LG2=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l2).*laguerre(p, l2, 2*r.^2/W^2).*exp(Gouy);



LG=abs((LG1+LG2).^2);
LG_n=LG./max(LG(:));
LG_a=angle(LG1);
imagesc(LG_n);

G=[1 1 1];
for i=0:50
    G=[G;[1-0.02*i 1 1-0.02*i]];
end
colormap(hot)
shading interp; lighting phong; view(2); axis equal; axis tight; axis off;

function y=laguerre(p,l,x)
y=zeros(p+1,1);
L=abs(l);
if p==0
    y=1;
else
for m=0:p
    y(p+1-m)=((-1).^m.*(factorial(p+L)))./(factorial(p-m).*factorial(L+m).*factorial(m));
end
end
y=polyval(y,x);
end
%function y=LGMode(p,l,rad)
%C=(2*factorial(p)./(pi*factorial(p+abs(l))))^0.5;
%Gouy=-1j*((-(l.*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k.*(rad.^2./(2*R))+k*z);
%y=C/W.*(2^0.5.*rad./W)^abs(l).*laguerre(p, l, 2.*rad^2./W^2).*exp((-1).*rad.^2./W^2).*exp(Gouy);
%end



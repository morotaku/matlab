clear all
close all
clc
l1=26;
l2=-26;
l3=24;
l4=-24;
l5=20;
l6=-20;
p1=0;
p2=0;
p3=0;
p4=0;
p5=0;
p6=0;
w=1;% micrrons
lam=640*10^-3; % micrrons
k=2*pi/lam;
Ly=10; % micrrons
Lx=10; % micrrons
N=1000;
dy=Ly/N;
M=1000;
dx=Lx/M;
z=2.0;
y=-Ly/2:dy:Ly/2-dy;
x=-Lx/2:dx:Lx/2-dx;
[x,y]=meshgrid(x,y);
%r=sqrt(x.^2+y.^2);
%[phi,r] = cart2pol(x,y);
[phi,r] = cart2pol(x,y);
zr=w^2*k/2;
Z=z/zr
R=z+zr^2/z;              
W=w*(1+(z/zr)^2)^0.5;  
%%
%%%%%%%%% LG1 %%%%%%%%%%%%%
C1=(2*factorial(p1)/(pi*factorial(p1+abs(l1))))^0.5;
Gouy1=1i*k*(z+r.^2/2/R)-1i*(2*p1+abs(l1)+1)*atan(z/zr);
LGpl1=C1/W*exp(1i*l1*phi).*(2^0.5/W*r).^abs(l1).*laguerre(p1,abs(l1),2/W^2*r.^2).*exp(-r.^2/W^2).*exp(Gouy1);
%%%%%%%%% LG2 %%%%%%%%%%%%%
C2=(2*factorial(p2)/(pi*factorial(p2+abs(l2))))^0.5;
Gouy2=1i*k*(z+r.^2/2/R)-1i*(2*p2+abs(l2)+1)*atan(z/zr);
LGpl2=C2/W*exp(1i*l2*phi).*(2^0.5/W*r).^abs(l2).*laguerre(p2,abs(l2),2/W^2*r.^2).*exp(-r.^2/W^2).*exp(Gouy2);
%%%%%%%%% LG3 %%%%%%%%%%%%%
C3=(2*factorial(p3)/(pi*factorial(p3+abs(l3))))^0.5;
Gouy3=1i*k*(z+r.^2/2/R)-1i*(2*p3+abs(l3)+1)*atan(z/zr);
LGpl3=C3/W*exp(1i*l3*phi).*(2^0.5/W*r).^abs(l3).*laguerre(p3,abs(l3),2/W^2*r.^2).*exp(-r.^2/W^2).*exp(Gouy3);
%%%%%%%%% LG4 %%%%%%%%%%%%%
C4=(2*factorial(p4)/(pi*factorial(p4+abs(l4))))^0.5;
Gouy4=1i*k*(z+r.^2/2/R)-1i*(2*p4+abs(l4)+1)*atan(z/zr);
LGpl4=C4/W*exp(1i*l4*phi).*(2^0.5/W*r).^abs(l4).*laguerre(p4,abs(l4),2/W^2*r.^2).*exp(-r.^2/W^2).*exp(Gouy4);
%%%%%%%%% LG5 %%%%%%%%%%%%%
C5=(2*factorial(p5)/(pi*factorial(p5+abs(l5))))^0.5;
Gouy5=1i*k*(z+r.^2/2/R)-1i*(2*p5+abs(l5)+1)*atan(z/zr);
LGpl5=C5/W*exp(1i*l5*phi).*(2^0.5/W*r).^abs(l5).*laguerre(p5,abs(l5),2/W^2*r.^2).*exp(-r.^2/W^2).*exp(Gouy5);
%%%%%%%%% LG6 %%%%%%%%%%%%%
C6=(2*factorial(p6)/(pi*factorial(p6+abs(l6))))^0.5;
Gouy6=1i*k*(z+r.^2/2/R)-1i*(2*p6+abs(l6)+1)*atan(z/zr);
LGpl6=C6/W*exp(1i*l6*phi).*(2^0.5/W*r).^abs(l6).*laguerre(p6,abs(l6),2/W^2*r.^2).*exp(-r.^2/W^2).*exp(Gouy6);
%%
LGpl=LGpl1+LGpl2+LGpl3+LGpl4+LGpl5+LGpl6;
LG=abs(LGpl).^2;
LGa=angle(LGpl);
imagesc(LGa);
colormap(hot)
shading interp; lighting phong; view(2); axis equal; axis tight; axis off;

function y=laguerre(n,l,x)
%LAGUERRE Laguerre polynomial.
%   Y=LAGUERRE(N,L,X) returns the generalized Laguerre polynomial of order 
%   N and index L of X, where N and L are a non-negative integers and is X 
%   a real scalar.
%
%   EXAMPLE:
%   y=laguerre(1,2,3.3)
%
%   See also templ, hermite

%   G.R.B.E. Römer 13-dec-1994
%   Revised G.R.B.E. Römer 20-sep-2010
%   Copyright 2008-2010 G.R.B.E. Römer, University of Twente

if prod(size(n))>1
    error('Laser Toolbox: HERMITE: N should be a non-negative scalar integer.');
end

if (n<0 || abs(round(n)-n)>0)
    abs(round(n)-n)>0
   error('Laser Toolbox: HERMITE: N should be a non-negative integer.');
end

if prod(size(l))>1
    error('Laser Toolbox: HERMITE: L should be a non-negative scalar integer.');
end

if (l<0 || abs(round(l)-l)>0)
    abs(round(n)-n)>0
   error('Laser Toolbox: HERMITE: L should be a non-negative integer.');
end

if n==0
   y=1;
elseif n==1,
   y=-1*x+1+l;
elseif n>1
   y= (2*n+l-1-x)./n.*laguerre(n-1,l,x) - (n+l-1)/n*laguerre(n-2,l,x);
end
end

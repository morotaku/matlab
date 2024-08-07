clear all
close all
clc
format shortEng
z=10
%%% x-y　coordinate
N=200;
L=50; 
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
%%% Topological charge
l=3;
%%% Radial index
n=3;


%%% Parameters(µm)
d=1000; %Incident beam size
lam=0.64; %Wave length
k=2*pi/lam; %Wave number
f=10^5; %Focus length
w0=1000;%Beam waist
w=lam*f/(pi*d)
zr=w^2*k/2 %Rayleigh length
s=0; %Distance from the input plane to lens

G0=1;
i=1j;

%ABCD matrix
A=-z/f
B=-(z/f)*s+f+z;
C=-1/f;
D=1-s/f;
R=z+zr^2./z; %Beam curvature
q0=1j*k*w0^2/2;
Q=(A+B/q0)/(C+D/q0);
W=w0.*sqrt(A.^2+(B.*lam./(pi*w0^2)).^2);
a=((k.*r./(2*B)).^2);
b=(1./w0.^2);
c=(1j*k.*A./(2.*B));
lag=laguerre(3,a./(b+c));
function y=laguerre(n,x)
y=0;
    if n==0
        y=1;
    else
        for m=0:n
            y=y+(factorial(n)^2*(-1)^m.*x.^m)./(factorial(m)*factorial(m)*factorial(n-m));
        end
    end
end

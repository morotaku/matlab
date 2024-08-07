clear all
close all
clc
%format shortEng
%%% Topological charge
l=0;
%%% Radial index
n=10;

%%% Parameters(µm)
w=1000; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
f=5000; %Focus length
s=0; %Distance from the input plane to lens
z=input('position: '); %Beam position

G0=1;
i=1j;

%ABCD matrix
A=1-z/f;
B=s+z-z*s/f;
%A=-z/f;
%B=-(z/f)*s+f+z;
C=-1/f;
D=1-s/f;
A=1;
B=z;
C=0;
D=1;

q0=1j*k*w^2/2;
Q=(A+B/q0)/(C+D/q0);
W=w*sqrt(A^2+(B*lam./(pi*w^2))^2);
%%%%%%%%%%%%%%%%Intinsity plot%%%%%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
L=5000; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);

%%%Mode plot
C1=G0*factorial(n)/2^n;
C2=0;
%A=(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1))*exp(-1j*k*z).*laguerre(m,2/w^2.*r.^2)%.*exp(-1j*k/(2*Q).*r.^2)
for m=0:n
    C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1))*exp(-1j*k*z).*laguerre(m,2/w^2.*r.^2).*exp(-1j*k/(2*Q).*r.^2);
    %C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1)).*laguerre(m,2/w^2.*r.^2).*exp(-1j*k/(2*Q).*r.^2);
end

lag=C1.*C2;
HGM=abs(lag.^2);
HGM_n=HGM/max(HGM(:));
imagesc(X,Y,HGM_n);

colormap("hot")
shading interp; lighting phong; view(2); axis equal; axis tight;% axis off;
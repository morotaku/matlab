clear all
close all
clc
format shortEng
%%% Topological charge
l=0;
%%% Radial index
n=3;

%%% Parameters(µm)
w=1000; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
f=10^5; %Focus length
%s=10^6; %Distance from the input plane to lens

G0=1;
i=1j;

z_inp=input('position: ');
z=z_inp;%+10^5;

%ABCD matrix
A=1-z/f;
B=z;
C=-1/f;
D=1;
%W=w*sqrt(A^2+(B*lam/(pi*w^2))^2);
%w2=w/(sqrt(1-(s/f)^2+(zr/f)^2))
w2=(lam*f/(pi*w))
zr2=pi*w2^2/lam
q0=1j*k*w2^2/2;
Q=(A+B/q0)/(C+D/q0);
%W=w2*(1+(z/zr2)^2)^0.5;
W=w2*sqrt(A^2+(B*lam/(pi*w2^2))^2)


%%%%%%%%%%%%%%%%Intinsity plot%%%%%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
L=1*10^2; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);

%%%Mode plot
C1=G0*factorial(n)/2^n;
C2=0;
for m=0:n
    C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1))*exp(-1j*k*z).*laguerre(m,(2.*r.^2)./W^2).*exp(-1j*k/(2*Q).*r.^2);
end

hgb=C1.*C2;
HGM=abs(hgb).^2;
HGM_n=HGM/max(HGM(:));
imagesc(X,Y,HGM_n);

colormap("hot")
shading interp; lighting phong; view(2); axis equal; axis tight;% axis off;

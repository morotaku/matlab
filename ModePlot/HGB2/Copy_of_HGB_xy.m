clear all
close all
clc
format shortEng
z=input('z=')
%%% x-y　coordinate
N=200;
L=20; 
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
%%% Beam index
n=10;
p=0;
l=1;


%%% Parameters(µm)
w0=10^5; %Incident beam size
lam=1.06; %Wave length
k=2*pi/lam; %Wave number
f=2*10^5; %Focus length
w=lam*f/(pi*w0) %Beam waist
zr=w^2*k/2 %Rayleigh length
s=2*10^6; %Distance from the input plane to lens

G0=1;
i=1j;

%ABCD matrix
A=-z/f;
B=-(z/f)*s+f+z;
C=-1/f;
D=1-s/f;
R=z+zr^2./z; %Beam curvature
q0=1j*k*w0^2/2;
Q=(A+B./q0)./(C+D./q0);
W=w*(1+(z/zr)^2)^0.5; %Beam size
%%%HGB
hgb=LGmode(p,l,r,phi,z,20.05,lam);
hgb_n=abs(hgb)./max(max(abs(hgb)));
HGB=hgb_n.^2;
HGB_n=HGB./max(HGB(:));

figure(1)
imagesc(X,Y,HGB_n,[0,1])
colormap("hot")
axis equal; axis tight;
figure(2)
plot(X,HGB_n(:,N/2))
%xlim([0 4000]);
%shading interp; lighting phong; view(2); 
%axis off;


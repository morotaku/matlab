clear all
close all
clc
format shortEng
%%% Topological charge
l=0;
%%% Radial index
p=10;

%%% Parameters(µm)
w=1; %Beam waist
lam=1.06; %Wave length

k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
z=input('position: '); %Beam position
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size


f=20000; %Focus length
s=100; %Distance from the input plane to lens
G0=1;
i=1j;

%ABCD matrix
A=-z/f;
B=(-z/f)*s+f+z;
C=-1/f;
D=1-s/f;


%%%%%%%%%%%%%%%%Intinsity plot%%%%%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N1=1001;
L=10000; %Display lange
X1=linspace(-L,L,N1);
Y1=linspace(-L,L,N1);
[x1,y1]=meshgrid(X1,Y1);
[phi2,r1] = cart2pol(x1,y1);

%%%Mode plot
C1=(i*A*k*G0*factorial(p)/(2*B*w^(2*p)))*((1/(w^2)+(i*k*A/(2*B)))^(-1-p))*exp(-i*k*z);
C2=exp(((-i*k*D.*(r1.^2))./(2*B)));
C3=exp(-((k.*r1)./(2*B)).^2./(1/w^2+(i*k*A/(2*B))));
lag=laguerre(p,(k.*r1./(2*B)).^2./((1/w^2)+(i*k*A/(2*B))));
hgb=C1.*C2.*C3.*lag;
HGB=abs(hgb).^2;
HGBn=HGB./(max(max(HGB)));
imagesc(X1,Y1,HGBn);
%im.max(0.8);


%%%%%%%%%%%%%%%%Polarization plot%%%%%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=30;
L=100; %Display lange
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);
[phi2,r2] = cart2pol(x2,y2);



colormap("hot")
shading interp; lighting phong; view(2); axis equal; axis tight;% axis off;
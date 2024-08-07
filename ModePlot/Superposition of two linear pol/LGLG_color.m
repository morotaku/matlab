clear all
close all
clc

%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=0.75;

%%% Topological charge
l1=1; %North pole
l2=-1; %South pole
%%% Radial index
p1=0;
p2=0;
%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude

w=1; %Beam waist
lam=0.64; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size


%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1000;
L=8; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);

[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
LG1_i=lgmode(p1,l1,r,phi1,z,w,lam);
LG2_i=lgmode(p2,l2,r,phi1,z,w,lam);

i_r=0;
for t = 0:20
    E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi1)).*exp(1j*t*pi/20).*LG1_i;  %ℓ=1
    E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi1)).*exp(1j*t*pi/20).*LG2_i; %ℓ=-1
    ex1_i=E1_i;
    ey1_i=0;
    ex2_i=0;
    ey2_i=E2_i;
    ex_i=ex1_i+ex2_i;
    ey_i=ey1_i+ey2_i;
    [Ex_i,Ey_i]=polarizer(angle,ex_i,ey_i);
    i_r=i_r+real(Ex_i.*conj(Ex_i)+Ey_i.*conj(Ey_i));
end
I=i_r./max(max(i_r));
imagesc([-L L],[-L L],I);
colormap('gray')
hold on
axis equal; axis off;
%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=11;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
lg1=lgmode(p1,l1,r2,phi2,z,w,lam);
lg2=lgmode(p2,l2,r2,phi2,z,w,lam);
phase_d=(l1.*phi2)-(l2.*phi2);
aa=sin(phase_d);




cnt1=0;
cnt2=0;
cnt3=0;
g1=[];
g2=[];
g3=[];

for i=1:N2^2
    if abs(aa(i))<10^(-10)
        cnt1=cnt1+1;
        g1(cnt1)=i;
        LG1_1(cnt1)=lg1(i);
        LG2_1(cnt1)=lg2(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);
    elseif aa(i)>0
        cnt2=cnt2+1;
        g2(cnt2)=i;
        LG1_2(cnt2)=lg1(i);
        LG2_2(cnt2)=lg2(i);
        x_2(cnt2)=x2(i);
        y_2(cnt2)=y2(i);
    elseif aa(i)<0
        cnt3=cnt3+1;
        g3(cnt3)=i;
        LG1_3(cnt3)=lg1(i);
        LG2_3(cnt3)=lg2(i);
        x_3(cnt3)=x2(i);
        y_3(cnt3)=y2(i);
    end
end

if cnt1==0
    LG1_1=0;
    LG2_1=0;
    x_1=0;
    y_1=0;
elseif cnt2==0
    LG1_2=0;
    LG2_2=0;
    x_2=0;
    y_2=0;
elseif cnt3==0
    LG1_3=0;
    LG2_3=0;
    x_3=0;
    y_3=0;
end

%%%Linear
for t = 0:120
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g1))).*exp(1j*t*pi/60).*LG1_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g1))).*exp(1j*t*pi/60).*LG2_1; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=E1_1;
    ey1_1=0;
    %%%|RHC> ℓ=-1
    ex2_1=0;
    ey2_1=E2_1;
    
    ex_1=real(ex1_1+ex2_1);
    ey_1=real(ey1_1+ey2_1);
    [Ex_1,Ey_1]=polarizer(angle,ex_1,ey_1);
    Ex1=Ex_1.*scale+x_1;
    Ey1=Ey_1.*scale+y_1;
    axis equal; %axis off;
    
    scatter(Ex1,Ey1,width,'.','green')
    pause(0.005)
end
hold on

%%%Left circular
for t = 0:120
    %%%North pole of Poincare sphere
    E1_2=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g2))).*exp(1j*t*pi/60).*LG1_2;
    %%%South pole
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g2))).*exp(1j*t*pi/60).*LG2_2; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_2=real(E1_2);
    ey1_2=0;
    %%%|RHC> ℓ=-1
    ex2_2=0;
    ey2_2=real(E2_2);
    
    ex_2=ex1_2+ex2_2;
    ey_2=ey1_2+ey2_2;
    [Ex_2,Ey_2]=polarizer(angle,ex_2,ey_2);
    Ex2=Ex_2.*scale+x_2;
    Ey2=Ey_2.*scale+y_2;
    axis equal;axis tight;% axis off;
    
    scatter(Ex2,Ey2,width,'.','red')
    pause(0.005)
end
hold on
fig=gcf;
fig.Color = 'none';
fig.InvertHardcopy = 'off';

%%%Right circular
for t = 0:120
    %%%North pole of Poincare sphere
    E1_3=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g3))).*exp(1j*t*pi/60).*LG1_3;
    %%%South pole
    E2_3=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g3))).*exp(1j*t*pi/60).*LG2_3; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_3=real(E1_3);
    ey1_3=0;
    %%%|RHC> ℓ=-1
    ex2_3=0;
    ey2_3=real(E2_3);
    
    ex_3=ex1_3+ex2_3;
    ey_3=ey1_3+ey2_3;
    [Ex_3,Ey_3]=polarizer(angle,ex_3,ey_3);
    Ex3=Ex_3.*scale+x_3;
    Ey3=Ey_3.*scale+y_3;
    axis equal;axis tight;% axis off;
    
    scatter(Ex3,Ey3,width,'.','blue')
    pause(0.005)
end
shading interp; lighting phong; view(2); axis tight;
axis off;


function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l).*laguerre(p, l, 2*r.^2/W^2).*exp(Gouy);
LG=abs(LG1.^2);
LG_n=LG./max(LG(:));
y=LG_n;
end

function y=lgmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l).*laguerre(p, l, 2*r.^2/W^2).*exp(Gouy);
LG=abs(LG1.^2);
LG_n=LG./max(LG(:));
y=LG1./max(max(LG1));
end

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

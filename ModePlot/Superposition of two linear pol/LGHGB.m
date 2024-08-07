clear all
close all
clc

%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
L=35; %Display lange
scale=L/10*0.9;
%%% Topological charge
l1=1; %North pole
n=1; %South pole
%%% Radial index
p1=0;
p2=0;
%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude

w=20.05; %Beam waist
w_h=1000;
lam=0.64; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
f=10^5;
s=0;
G0=1;


%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1000;

X=linspace(-L,L,N);
Y=linspace(-L,L,N);

[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
LG1_i=LGmode(p1,l1,r,phi1,z,w,lam);
hgb1_i=f_hgb(n,f,z,lam,w_h,G0,r);

i_r=0;
for t = 0:20
    E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi1)).*exp(1j*t*pi/20).*LG1_i;  %ℓ=1
    E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/20).*hgb1_i; %ℓ=-1
    ex1_i=real(E1_i);
    ey1_i=0;
    ex2_i=0;
    ey2_i=real(E2_i);
    ex_i=ex1_i+ex2_i;
    ey_i=ey1_i+ey2_i;
    [Ex_i,Ey_i]=polarizer(angle,ex_i,ey_i);
    i_r=i_r+real(sqrt(Ex_i.^2+Ey_i.^2));
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
lg1=LGmode(p1,l1,r2,phi2,z,w,lam);
hgb1=f_hgb(n,f,z,lam,w_h,G0,r2);


%%%Linear
for t = 0:120
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2)).*exp(1j*t*pi/60).*lg1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*hgb1; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=real(E1_1);
    ey1_1=0;
    %%%|RHC> ℓ=-1
    ex2_1=0;
    ey2_1=real(E2_1);
    
    ex_1=ex1_1+ex2_1;
    ey_1=ey1_1+ey2_1;
    [Ex_1,Ey_1]=polarizer(angle,ex_1,ey_1);
    Ex1=Ex_1.*scale+x2;
    Ey1=Ey_1.*scale+y2;
    axis equal; %axis off;
    
    scatter(Ex1,Ey1,width,'.','red')
    pause(0.005)
end
shading interp; lighting phong; view(2); axis tight;
axis off;

function y=f_hgb(n,f,z,lam,w0,G0,r)
    %ABCD matrix
    s=0;
    A=-z./f;
    B=-(z./f).*s+f+z;
    C=-1/f;
    D=1-s/f;
    k=2*pi/lam; %Wave number
    q0=1j*k.*w0.^2./2;
    %Q=(A+B./q0)./(C+D./q0);
    %W=w0.*sqrt(A.^2+(B.*lam./(pi*w0^2)).^2);

    %%%HGB
    C1=1j*k*G0*factorial(n)./(2*w0^(2*n).*B).*((1./w0.^2)+(1j*k.*A./(2.*B))).^(-1-n).*exp(-1j*k.*z);
    C2=exp(-1j*k.*(D.*r.^2)./(2.*B));
    C3=exp(-(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B)))).*laguerre(n,(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B))));
    hgb=C1.*C2.*C3;
    hgb_n=abs(hgb)./max(max(abs(hgb)));
    HGB=hgb_n.^2;
    HGB_n=HGB./max(HGB(:));
    y=HGB_n;
end

function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l).*a_laguerre(p, l, 2*r.^2/W^2).*exp(Gouy);
LG=abs(LG1.^2);
LG_n=LG./max(LG(:));
y=LG_n;
end

function y=a_laguerre(p,l,x)
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
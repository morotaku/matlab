clear all
close all
clc

%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
L=450; %Display lange
scale=L/10*0.7;
%%% Topological charge
l1=4; %North pole
n=2; %South pole
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
N=1001;

X=linspace(-L,L,N);
Y=linspace(-L,L,N);

[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);
r=sqrt(x.^2+y.^2);
phi1=atan2(y,x);

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
N2=21;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);
r2=sqrt(x2.^2+y2.^2);
phi2=atan2(y2,x2);

%%% LG mode
lg1=LGmode(p1,l1,r2,phi2,z,w,lam);
hgb1=f_hgb(n,f,z,lam,w_h,G0,r2);
phase_d=(l1.*phi2);
aa=sin(phase_d);



cnt1=0;
cnt2=0;
cnt3=0;
g1=[];
g2=[];
g3=[];

for i=1:N2^2
    if abs(aa(i))<0.1%10^(-10)
        cnt1=cnt1+1;
        g1(cnt1)=i;
        LG1_1(cnt1)=lg1(i);
        HGB2_1(cnt1)=hgb1(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);
    elseif aa(i)>0
        cnt2=cnt2+1;
        g2(cnt2)=i;
        LG1_2(cnt2)=lg1(i);
        HGB2_2(cnt2)=hgb1(i);
        x_2(cnt2)=x2(i);
        y_2(cnt2)=y2(i);
    elseif aa(i)<0
        cnt3=cnt3+1;
        g3(cnt3)=i;
        LG1_3(cnt3)=lg1(i);
        HGB2_3(cnt3)=hgb1(i);
        x_3(cnt3)=x2(i);
        y_3(cnt3)=y2(i);
    end
end

if cnt1==0
    LG1_1=0;
    HGB2_1=0;
    x_1=0;
    y_1=0;
elseif cnt2==0
    LG1_2=0;
    HGB2_2=0;
    x_2=0;
    y_2=0;
elseif cnt3==0
    LG1_3=0;
    HGB2_3=0;
    x_3=0;
    y_3=0;
end

%%%Linear
for t = 0:120
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g1))).*exp(1j*t*pi/60).*LG1_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*HGB2_1; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=real(E1_1);
    ey1_1=0;
    %%%|RHC> ℓ=-1
    ex2_1=0;
    ey2_1=real(E2_1);
    
    ex_1=ex1_1+ex2_1;
    ey_1=ey1_1+ey2_1;
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
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*HGB2_2; %ℓ=-1
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
    E2_3=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*HGB2_3; %ℓ=-1
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
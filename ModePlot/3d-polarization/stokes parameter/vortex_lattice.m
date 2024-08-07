clear all
close all
clc
L=2.75; %Display lange
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=L/10.*1.2;
%%% Polarizer angle 
ang=input('polarizer angle: '); %thita=-1 → no polarizer

%%% Topological charge
n1=0; %North pole
n2=2; %South pole
%%% Radial index
m1=1;
m2=0;

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
N=1001;

X=linspace(-L,L,N);
Y=linspace(-L,L,N);

[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
hg1=HGmode(n1,m1,x,y,w,r,z,R,k);
hg2=HGmode(n2,m2,x,y,w,r,z,R,k);
LG1_i=hg1;
LG2_i=hg2;

i_r=0;
i_x=0;
i_y=0;
for t = 0:40
    ex=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/20).*LG1_i;%.*1.5;  %ℓ=1
    ey=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/20).*LG2_i.*1j; %ℓ=-1
    [Ex_i,Ey_i]=polarizer(ang,ex,ey);
    %i_r=i_r+real(sqrt(Ex_i.^2+Ey_i.^2));
    i_x=i_x+Ex_i;
    i_y=i_y+Ey_i;
end
U=Iy;
U_n=U./max(U(:));
I=i_r./max(max(i_r));

figure(1)
imagesc([-L L],[-L L],U_n);
colormap('gray')
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
hg1=HGmode(n1,m1,x,y,w,r,z,R,k);
hg2=HGmode(n2,m2,x,y,w,r,z,R,k);
hg1=hgmode(0,4,x2,y2,w,r2,z,R,k);
hg2=HGmode(n2,m2,x2,y2,w,r2,z,R,k);
LG_q=hg1./max(max(abs(hg1)));
HGB_q=hg2./max(max(abs(hg2)));

cnt1=0;
cnt2=0;
g1=[];
g2=[];

for i=1:N2^2
    if abs(LG_q(i))>=abs(HGB_q(i))
        cnt1=cnt1+1;
        g1(cnt1)=i;
        LG_1(cnt1)=LG_q(i);
        HGB_1(cnt1)=HGB_q(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);

        
    else
        cnt2=cnt2+1;
        g2(cnt2)=i;
        LG_2(cnt2)=LG_q(i);
        HGB_2(cnt2)=HGB_q(i);
        x_2(cnt2)=x2(i);
        y_2(cnt2)=y2(i);
    end
end

if cnt1==0
    LG_1=0;
    HGB_1=0;
    x_1=0;
    y_1=0;
elseif cnt2==0
    LG_2=0;
    HGB_2=0;
    x_2=0;
    y_2=0;
end



%%%Right circular
for t = 0:120
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/60).*LG_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*HGB_1.*exp(1j*pi/2); %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=real(E1_1);
    ey1_1=0;
    %%%|RHC> ℓ=-1
    ex2_1=0;
    ey2_1=real(E2_1);
    
    ex_1=ex1_1+ex2_1;
    ey_1=ey1_1+ey2_1;
    [Ex_1,Ey_1]=polarizer(ang,ex_1,ey_1);
    Ex1=Ex_1.*scale+x_1;
    Ey1=Ey_1.*scale+y_1;
    axis equal; %axis off;
    hold on
    scatter(Ex1,Ey1,width,'.','red')
    pause(0.005)
end
hold on
fig=gcf;
fig.Color = 'none';
fig.InvertHardcopy = 'off';

%%%Left circular
for t = 0:120
    %%%North pole of Poincare sphere
    E1_2=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/60).*LG_2;
    %%%South pole
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*HGB_2.*exp(1j*pi/2); %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_2=real(E1_2);
    ey1_2=0;
    %%%|RHC> ℓ=-1
    ex2_2=0;
    ey2_2=real(E2_2);
    
    ex_2=ex1_2+ex2_2;
    ey_2=ey1_2+ey2_2;
    [Ex_2,Ey_2]=polarizer(ang,ex_2,ey_2);
    Ex2=Ex_2.*scale+x_2;
    Ey2=Ey_2.*scale+y_2;
    axis equal;axis tight;% axis off;
    
    scatter(Ex2,Ey2,width,'.','blue')
    pause(0.005)
end
shading interp; lighting phong; view(2); axis tight;
axis off;



function y=HGmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-(r.^2./W^2)-(1j*k.*r.^2/(2*R)));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u./max(max(abs(u)));%./max(max(abs(u)));
end
function y=hgmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-(r.^2./W^2)-(1j*k.*r.^2/(2*R)));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));
end

function y=HermitePol(n,x)
Hn=0;
for m=0:floor(n/2)
    Hn=Hn+((-1)^m.*(2.*x).^(n-2*m))./(factorial(m)*factorial(n-2*m));
end
Hn=Hn.*factorial(n);
y=Hn;
end

function [x,y]=polarizer(thita,ex,ey)
if thita==-1
    x=ex;
    y=ey;
else
    x=cos(thita)^2.*ex+cos(thita)*sin(thita).*ey;
    y=sin(thita)*cos(thita).*ex+sin(thita)^2.*ey;
end
end
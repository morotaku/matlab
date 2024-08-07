clear all
close all
clc
i=1j;
L=9;
pass='C:\Users\maila\OneDrive\デスクトップ\保存先\4\270.emf';
%%%Figure parameter
width=10;
scale=L/2;
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%

%%% Topological charge
l1=1; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
p2=0;

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude強度比

%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of LG mode
w1=5; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2 %Rayleigh length



%%% parameter of HGB
w2=5; %Incident Beam size
lam=0.63; %Wave length
zr2=w2^2*k/2; %Rayleigh length
f=10^5; %Focus length
s=zr2; %Distance from the input plane to lens

w2_2=w2/(sqrt(1-(s/f)^2+(zr1/f)^2)); %Beam waist after focus
zr2_2=pi*w2_2^2/lam; %Rayleigh length after focus

z1=input('position: '); %Beam position
z2=z1;


G0=1;


%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
%L=35; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
LG_i1=LGmode(p1,l1,r,phi1,z1,w1,lam);

%%% HGB
LG_i2=LGmode(p2,l2,r,phi1,z1,w2,lam);

i_x=0; i_y=0;
i_r=0;
for t = 0:40
    E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/20).*LG_i1;  %ℓ=1
    E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/20).*LG_i2; %ℓ=-1
    ex1_i=E1_i;
    ey1_i=-1j.*E1_i;
    ex2_i=E2_i;
    ey2_i=i.*E2_i;
    ex_i=ex1_i+ex2_i;
    ey_i=ey1_i+ey2_i;
    [Ex_i,Ey_i]=polarizer(angle,ex_i,ey_i);
    %i_r=i_r+real((Ex_i.*conj(Ex_i)+Ey_i.*conj(Ey_i)));
    i_x=i_x+Ex_i;
    i_y=i_y+Ey_i;
end
I=real(i_x.*conj(i_x))+real(i_y.*conj(i_y));
I_n=I./max(max(I));
%I=i_r./max(max(i_r));

%axis off;
figure(1)
imagesc(X,Y,I_n);
xticks(-L:L/5:L)
yticks(-L:L/5:L)
colormap('gray')
axis equal; axis off;
shading interp; lighting phong; view(2); axis tight;
hold on
%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=13;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
LG_n1=LGmode(p1,l1,r2,phi2,z1,w1,lam).*exp(1j.*pi/2);

%%% HGB
LG_n2=LGmode(p2,l2,r2,phi2,z1,w2,lam);



cnt1=0;
cnt2=0;
g1=[];
g2=[];

for i=1:N2^2
    if cos(thita/2).*abs(LG_n1(i))>=sin(thita/2).*abs(LG_n2(i))
        cnt1=cnt1+1;
        g1(cnt1)=i;
        LG_r1(cnt1)=LG_n1(i);
        LG_l1(cnt1)=LG_n2(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);

        
    else
        cnt2=cnt2+1;
        g2(cnt2)=i;
        LG_r2(cnt2)=LG_n1(i);
        LG_l2(cnt2)=LG_n2(i);
        x_2(cnt2)=x2(i);
        y_2(cnt2)=y2(i);
    end
end



%%%Right circular
for t = 0:120
    %%%North pole of Poincare sphere
    Er_1=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/60).*LG_r1;  %ℓ=1
    %%%South pole
    El_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*LG_l1; %ℓ=-1
    %%%|LHC> ℓ=1
    exr_1=real(Er_1);
    eyr_1=real(1j.*Er_1);
    %%%|RHC> ℓ=-1
    exl_1=real(El_1);
    eyl_1=real(-1j.*El_1);
    
    ex_1=exr_1+exl_1;
    ey_1=eyr_1+eyl_1;
    [Ex_1,Ey_1]=polarizer(angle,ex_1,ey_1);
    Ex1=Ex_1.*scale+x_1;
    Ey1=Ey_1.*scale+y_1;
    axis equal; axis off;
    
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
    Er_2=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/60).*LG_r2;%.*exp(1j.*(l1.*phi2(g2)))
    %%%South pole
    El_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*LG_l2; %.*exp(1j.*(l2.*phi2(g2)))
    %%%|LHC> ℓ=1
    ex1_2=real(Er_2);
    ey1_2=real(1j.*Er_2);
    %%%|RHC> ℓ=-1
    ex2_2=real(El_2);
    ey2_2=real(-1j.*El_2);
    
    ex_2=ex1_2+ex2_2;
    ey_2=ey1_2+ey2_2;
    [Ex_2,Ey_2]=polarizer(angle,ex_2,ey_2);
    Ex2=Ex_2.*scale+x_2;
    Ey2=Ey_2.*scale+y_2;
    axis equal;axis tight; axis off;
    
    scatter(Ex2,Ey2,width,'.','blue')
    pause(0.005)
end


%ax.Color = 'none';
%ax.Box = 'off';
%ax.XColor = 'none';
%ax.YColor = 'none';
%saveas(fig, pass);

function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2./z; %Beam curvature
W=w.*(1+(z./zr).^2).^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*(r.^2./(2*R))+k.*z);
lg=C./W.*exp((-1).*r.^2./W.^2).*(2^0.5*r./W).^abs(l).*laguerre(p, l, 2*r.^2./W.^2).*exp(Gouy);
LG=abs(lg.^2);
%y=LG./max(LG(:));
y=lg;
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

function [x,y]=polarizer(thita,ex,ey)
if thita==-1
    x=ex;
    y=ey;
else
    x=cos(thita)^2.*ex+cos(thita)*sin(thita).*ey;
    y=sin(thita)*cos(thita).*ex+sin(thita)^2.*ey;
end
end


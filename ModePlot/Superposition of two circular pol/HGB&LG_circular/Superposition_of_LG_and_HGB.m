clear all
close all
clc
i=1j;
L=75;
pass='C:\Users\maila\OneDrive\デスクトップ\保存先\4\270.emf';
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%

%%% Topological charge
l1=1; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
n=1;

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude強度比

%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of LG mode
w1=20.5; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2 %Rayleigh length



%%% parameter of HGB
w2=1000; %Incident Beam size
lam=0.63; %Wave length
zr2=w2^2*k/2; %Rayleigh length
f=100000; %Focus length
s=zr2; %Distance from the input plane to lens

w2_2=w2/(sqrt(1-(s/f)^2+(zr1/f)^2)); %Beam waist after focus
zr2_2=pi*w2_2^2/lam; %Rayleigh length after focus

z1=input('position: '); %Beam position
z2=z1+f;


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
LG_n=LGmode(p1,l1,r,phi1,z1,w1,lam);

%%% HGB
HGB=f_hgb(n,f,z2,lam,w2,G0,r);

i_r=0;
for t = 0:20
    E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi1)).*exp(1j*t*pi/20).*LG_n;  %ℓ=1
    E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi1)).*exp(1j*t*pi/20).*HGB; %ℓ=-1
    ex1_i=real(E1_i);
    ey1_i=real(-1j.*E1_i);
    ex2_i=real(E2_i);
    ey2_i=real(1j.*E2_i);
    ex_i=ex1_i+ex2_i;
    ey_i=ey1_i+ey2_i;
    [Ex_i,Ey_i]=polarizer(angle,ex_i,ey_i);
    i_r=i_r+real(sqrt(Ex_i.^2+Ey_i.^2));
end

I=i_r./max(max(i_r));
imagesc(X,Y,I);
xticks(-L:L/5:L)
yticks(-L:L/5:L)
colormap('gray')
hold on
axis equal; %axis off;
%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=11;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
LG_q=LGmode(p1,l1,r2,phi2,z1,w1,lam);

%%% HGB
HGB_q=f_hgb(n,f,z2,lam,w2,G0,r2);



cnt1=0;
cnt2=0;
g1=[];
g2=[];

for i=1:N2^2
    if sin(thita/2)*abs(LG_q(i))>=cos(thita/2)*abs(HGB_q(i))
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

L2=11;


%%%Figure parameter
width=10;
scale=L/10*0.9;

%%%Right circular
for t = 0:120
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g1))).*exp(1j*t*pi/60).*LG_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g1))).*exp(1j*t*pi/60).*HGB_1; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=real(E1_1);
    ey1_1=real(-1j.*E1_1);
    %%%|RHC> ℓ=-1
    ex2_1=real(E2_1);
    ey2_1=real(1j.*E2_1);
    
    ex_1=ex1_1+ex2_1;
    ey_1=ey1_1+ey2_1;
    [Ex_1,Ey_1]=polarizer(angle,ex_1,ey_1);
    Ex1=Ex_1.*scale+x_1;
    Ey1=Ey_1.*scale+y_1;
    axis equal; %axis off;
    
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
    E1_2=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g2))).*exp(1j*t*pi/60).*LG_2;
    %%%South pole
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g2))).*exp(1j*t*pi/60).*HGB_2; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_2=real(E1_2);
    ey1_2=real(-1j.*E1_2);
    %%%|RHC> ℓ=-1
    ex2_2=real(E2_2);
    ey2_2=real(1j.*E2_2);
    
    ex_2=ex1_2+ex2_2;
    ey_2=ey1_2+ey2_2;
    [Ex_2,Ey_2]=polarizer(angle,ex_2,ey_2);
    Ex2=Ex_2.*scale+x_2;
    Ey2=Ey_2.*scale+y_2;
    axis equal;axis tight;% axis off;
    
    scatter(Ex2,Ey2,width,'.','blue')
    pause(0.005)
end


%ax.Color = 'none';
%ax.Box = 'off';
%ax.XColor = 'none';
%ax.YColor = 'none';
%saveas(fig, pass);
clear all
close all
clc
L=3; %Display lange
p=3;m=3;ec=3;
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=L/10;
%%% Polarizer angle 
ang=input('polarizer angle(thita=-1 → no polarizer): '); 

%%% Topological charge
n1=0; %North pole
n2=3; %South pole
%%% Radial index
m1=3;
m2=0;


%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude

w=1; %Beam waist
lam=1; %Wave length
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
hg1=hgmode_45(4,x,y,w,r,z,R,k);
HG1=abs(hg1.^2);
HG1_n=HG1/max(max(HG1));
%hg2=hgmode(n2,m2,x,y,w,r,z,R,k);
[hg2,Xin,Yin]=Ince_Gaussian(15e-3,1001,0,p,m,ec,5e-3,(2*pi/lam),0);
HG2=abs(hg2.*conj(hg2));
HG2=rot90(HG2,-1);
HG2_n=HG2/max(max(HG2));

i_r=0;

i_x=0;
i_y=0;
for t = 0:40
    ex=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/20).*hg1;%.*1.5;  %ℓ=1
    ey=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/20).*hg2.*1j; %ℓ=-1
    [Ex_i,Ey_i]=polarizer(ang,ex,ey);
    i_r=i_r+real(ex.*conj(ex)+ey.*conj(ey));

%     i_x=i_x+Ex_i;
%     i_y=i_y+Ey_i;
%     i_x=i_x+Ex_i;
%     i_y=i_y+Ey_i;
end
%U=real((i_x).*conj(i_x))+real((i_y).*conj(i_y));
U=i_r;
U_n=U./max(U(:));
U_n=HG1_n+HG2_n;

imagesc([-L L],[-L L],U_n);
%imagesc([-L L],[-L L],II);
colormap('gray')
axis equal; %axis off;
hold on
%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=17;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
hg1=HGmode(n1,m1,x2,y2,w,r2,z,R,k);
hg1=hgmode_45(4,x2,y2,w,r2,z,R,k);
HG1=abs(hg1.^2);
lg1=HG1/max(max(HG1));
hg2=HGmode(n2,m2,x2,y2,w,r2,z,R,k);
[hg2,Xin,Yin]=Ince_Gaussian(15e-3,17,0,p,m,ec,5e-3,(2*pi/lam),0);
hg2=rot90(hg2./max(max(abs(hg2))));
HG2=abs(hg2.^2);
lg2=HG2/max(max(HG2));

aa=angle(hg1);
bb=angle(hg2.*exp(1j*pi/2));
cc=bb-aa;
dd=sin(cc);
% figure(2)
% imagesc(dd)

cnt1=0;
cnt2=0;
g1=[];
g2=[];

for i=1:N2^2
    if sin(cc(i))>0
        cnt1=cnt1+1;
        g1(cnt1)=i;
        LG_1(cnt1)=lg1(i);
        HGB_1(cnt1)=lg2(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);

        
    else
        cnt2=cnt2+1;
        g2(cnt2)=i;
        LG_2(cnt2)=lg1(i);
        HGB_2(cnt2)=lg2(i);
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

function y=HGmode_45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        if rem(i,2)==1
            hg=hg-nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        else
            hg=hg+nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        end
    end
    y=hg./(sqrt(2)^n);
    %y=hg./max(max(abs(hg)));
end

function y=HGmode45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        hg=hg+nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
    end
    y=hg./(sqrt(2)^n);
end

function y=hgmode_45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        if rem(i,2)==1
            hg=hg-nchoosek(n,i).*hgmode(i,n-i,x,y,w,r,z,R,k);
        else
            hg=hg+nchoosek(n,i).*hgmode(i,n-i,x,y,w,r,z,R,k);
        end
    end
    y=hg./(sqrt(2)^n);
    %y=hg./max(max(abs(hg)));
end

function y=hgmode45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        hg=hg+nchoosek(n,i).*hgmode(i,n-i,x,y,w,r,z,R,k);
    end
    y=hg./(sqrt(2)^n);
end

function y=HGmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-(r.^2./W^2)-(1j*k.*r.^2/(2*R)));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));
end

function y=hgmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-((x.^2+y.^2)./W^2));%-(1j*k.*(x.^2+y.^2)./(2*R))).*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));%./max(max(abs(u)));
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


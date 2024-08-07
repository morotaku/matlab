clear all
close all
clc

color='jet';
L=3;%1300;
pass='C:\Users\maila\OneDrive\デスクトップ\保存先\4\270.emf';
%%%Figure parameter
width=10;
scale=L/10*0.9;

%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%% Topological charge
l1=1; %North pole
l2=0; %South pole1
%%% Radial index
p1=0;
p2=0;

%%% Polarizer angle 
ang=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of LG mode
w1=1000; %Beam waist
w2=w1; 
lam=0.532; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2 %Rayleigh length

%z=input('position: '); %Beam position

%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N0=1001;
L0=1300;
%L=35; %Display lange
X0=linspace(-L0,L0,N0);
Y0=linspace(-L0,L0,N0);
[x0,y0]=meshgrid(X0,Y0);
z0=0;
%%% Polar coodinate
[phi0,r0] = cart2pol(x0,y0);
%%% LG mode
LG_i1=LGmode(p1,l1,r0,phi0,z0,w1,lam);

%%% HGB
LG_i2=LGmode(p2,l2,r0,phi0,z0,w2,lam);

%%%%%%%%%%%%%%Input beam plot%%%%%%%%%%%%%%%%%%%%
% Il=real(LG_i1.*conj(LG_i1)); %Input left circularly polarized field intensity
% Ir=real(LG_i2.*conj(LG_i2)); %Input right circularly polarized field intensity
%I=Il./max(max(Il))+Ir./max(max(Ir));%+Iz./max(max(Iz)); %Total input field intensity

ix_i=(LG_i1+LG_i2);
iy_i=1j.*(LG_i1-LG_i2);
Ix_i=real(ix_i.*conj(ix_i));
Iy_i=real(iy_i.*conj(iy_i));
I=Ix_i./max(max(Ix_i))+Iy_i./max(max(Iy_i));%+Iz./max(max(Iz)); %Total input field intensity


f1=figure(1);
imagesc(X0,Y0,I);
f1.Position(1:4) = [320 100 840 630];
xticks(-L0:L0/5:L0);yticks(-L0:L0/5:L0)
xlim([-L0 L0]);ylim([-L0 L0])
fontsize(1,15,"points")
axis tight;axis equal;
colormap(color)
colorbar
hold on

% ix=LG_i1+LG_i2;
% Ix=real(ix.*conj(ix));
% Ix_n=Ix./max(max(Ix));

% figure(2)
% imagesc(X0,Y0,I2);
% xticks(-L:L/5:L);yticks(-L:L/5:L)
% xlim([-1300 1300]);ylim([-1300 1300])
% axis tight;axis equal;
% colormap(color)
% colorbar
% 
% iy=1j.*LG_i1-1j.*LG_i2;
% Iy=real(iy.*conj(iy));
% Iy_n=Iy./max(max(Iy));

% figure(3)
% imagesc(X,Y,Iy_n);
% xticks(-L:L/5:L);yticks(-L:L/5:L)
% xlim([-1300 1300]);ylim([-1300 1300])
% axis tight;axis equal;
% colormap(color)
% colorbar


%%%%%%%%%%%%%%Focused beam plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
L=1.5;
%L=35; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi,r] = cart2pol(x,y);

%%% Parameter
NA=0.95;
f=5*10^3;
w0=0.61*lam/NA
zr=pi*w0^2/lam
z=zr.*0;

alpha=asin(NA);
thitaP=phi;


u=k*z*sin(alpha)^2;
v=k.*sqrt(x.^2+y.^2).*sin(alpha);

lx=LGmode(p2,l1,r,phi,0,w1,lam)+LGmode(p2,0,r,phi,0,w2,lam);
ly=1j.*LGmode(p2,l1,r,phi,0,w1,lam)-1j.*LGmode(p2,0,r,phi,0,w2,lam);
Ax=pi.*f.*lx./lam;
Ay=pi.*f.*ly./lam;

I0=0;I1=0;I2=0;
% for t=0:29
%     w=alpha/30;
%     thita=w*t;
%     thita2=w*(t+1);
%     a0=sqrt(cos(thita)).*sin(thita).*(1+cos(thita)).*besselj(0,v.*sin(thita)./sin(alpha)).*exp(1j.*u.*cos(thita)./sin(alpha).^2);
%     b0=sqrt(cos(thita2)).*sin(thita2).*(1+cos(thita2)).*besselj(0,v.*sin(thita2)./sin(alpha)).*exp(1j.*u.*cos(thita2)./sin(alpha).^2);
%     a1=sqrt(cos(thita)).*sin(thita).^2.*besselj(1,v.*sin(thita)./sin(alpha)).*exp(1j.*u.*cos(thita)./sin(alpha).^2);
%     b1=sqrt(cos(thita2)).*sin(thita2).^2.*besselj(1,v.*sin(thita2)./sin(alpha)).*exp(1j.*u.*cos(thita2)./sin(alpha).^2);
%     a2=sqrt(cos(thita)).*sin(thita).*(1-cos(thita)).*besselj(2,v.*sin(thita)./sin(alpha)).*exp(1j.*u.*cos(thita)./sin(alpha).^2);
%     b2=sqrt(cos(thita2)).*sin(thita2).*(1-cos(thita2)).*besselj(2,v.*sin(thita2)./sin(alpha)).*exp(1j.*u.*cos(thita2)./sin(alpha).^2);
%     I0=I0+(a0+b0).*w./2;
%     I1=I1+(a1+b1).*w./2;
%     I2=I2+(a2+b2).*w./2;
% end
for t=0:30
    w=alpha/30;
    thita=w*t;
    thita2=w*(t+1);
    a0=sqrt(cos(thita)).*sin(thita).*(1+cos(thita)).*besselj(0,v.*sin(thita)./sin(alpha)).*exp(1j.*u.*cos(thita)./sin(alpha).^2);
    a1=sqrt(cos(thita)).*sin(thita).^2.*besselj(1,v.*sin(thita)./sin(alpha)).*exp(1j.*u.*cos(thita)./sin(alpha).^2);
    a2=sqrt(cos(thita)).*sin(thita).*(1-cos(thita)).*besselj(2,v.*sin(thita)./sin(alpha)).*exp(1j.*u.*cos(thita)./sin(alpha).^2);
    I0=I0+(a0).*w;
    I1=I1+(a1).*w;
    I2=I2+(a2).*w;
end

I0_r=real(I0);
I1_r=real(I1);
I2_r=real(I2);

I0_i=imag(I0);
I1_i=imag(I1);
I2_i=imag(I2);

px=(Ax.*I0_i)+I2_i.*(Ax.*cos(2.*phi)+Ay.*sin(2.*phi));
py=(Ay.*I0_i)+I2_i.*(Ax.*sin(2.*phi)-Ay.*cos(2.*phi));
pz=-2.*I1_r.*(Ax.*cos(phi)+Ay.*sin(phi));

qx=-((Ax.*I0_r)+I2_r.*(Ax.*cos(2.*phi)+Ay.*sin(2.*phi)));
qy=-((Ay.*I0_r)+I2_r.*(Ax.*sin(2.*phi)-Ay.*cos(2.*phi)));
qz=-2.*I1_i.*(Ax.*cos(phi)+Ay.*sin(phi));

%%% intensity plot
i_x=real(px.*conj(px))+real(qx.*conj(qx));
i_y=real(py.*conj(py))+real(qy.*conj(qy));
i_z=real(pz.*conj(pz))+real(qz.*conj(qz));
i_t=i_x+i_y+i_z;

I_x=i_x./max(max(i_z));
I_y=i_y./max(max(i_z));
I_z=i_z./max(max(i_z));
I_t=I_x+I_y+I_z;

%%% phase plot
i_za=angle(pz+qz);
i_xa=angle(px+qx);
i_ya=angle(py+qy);


f2=figure(2);
hold on
imagesc(X,Y,I_x);
%imagesc(X,Y,i_xa);
f2.Position(1:4) = [320 100 840 630];
fontsize(2,15,"points")
xlim([-1000 1000]);ylim([-1300 1300])
xticks(-L:L/5:L);yticks(-L:L/5:L)
colormap(color)
colorbar
axis tight;axis equal;

f3=figure(3);
imagesc(X,Y,I_y);
%imagesc(X,Y,i_ya);
f3.Position(1:4) = [320 100 840 630];
fontsize(3,15,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal;

f4=figure(4);
imagesc(X,Y,I_z);
%imagesc(X,Y,i_za);
f4.Position(1:4) = [320 100 840 630];
fontsize(4,15,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal;

f5=figure(5);
imagesc(X,Y,i_t);
%imagesc(X,Y,i_xa+i_ya+i_za);
f5.Position(1:4) = [320 100 840 630];
fontsize(5,15,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal;




%%%%%%%%%%%%%%Function%%%%%%%%%%%%%%%%%%%%
function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2./z; %Beam curvature
W=w.*(1+(z./zr).^2).^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*(r.^2./(2*R))+k.*z);
LG1=C./W.*exp((-1).*r.^2./W.^2).*(2^0.5*r./W).^abs(l).*laguerre(p, l, 2*r.^2./W.^2).*exp(Gouy);
LG=abs(LG1.^2);
%y=LG./max(LG(:));
y=LG1./abs(max(max(LG1)));

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
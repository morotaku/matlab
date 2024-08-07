%%%Numerical study of spin–orbit interaction of light in nonparaxial focusing of Gaussian beams
clear all
close all
clc

color='jet';
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%% Topological charge
n1=0; %North pole
n2=1; %South pole1
%%% Radial index
m1=1;
m2=0;
%%% Input beam intensity and phase
ex_p=1;
ey_p=0;%1j

%%% Polarizer angle 
ang=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of Input beam
w1=1000;
w2=w1; 
lam=0.532; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2 %Rayleigh length
z=0;
R=z+zr1^2./z;
%z=input('position: '); %Beam position

%% Input beam image
%%% x-y　coordinate
N0=301;
L0=5000;
%L=35; %Display lange
X0=linspace(-L0,L0,N0);
Y0=linspace(-L0,L0,N0);
[x0,y0]=meshgrid(X0,Y0);
z0=0;
%%% Polar coodinate
[phi0,r0] = cart2pol(x0,y0);
%%% LG mode
%LG_i1=HGmode(n1,m1,x0,y0,w1,r0,0,R,k);
HG1=HGmode(n1,m1,x0,y0,w1,z,R,k);
%%% HGB
%LG_i2=1j.*HGmode(n2,m2,x0,y0,w2,r0,0,R,k);
HG2=HGmode(n2,m2,x0,y0,w1,z,R,k);

ix_i=HG1.*ex_p; 
iy_i=HG2.*ey_p;
Ix_i=real(ix_i.*conj(ix_i));
Iy_i=real(iy_i.*conj(iy_i));
I=Ix_i./max(max(Ix_i+0.0001))+Iy_i./max(max(Iy_i+0.0001));%+Iz./max(max(Iz)); %Total input field intensity


f1=figure(1);
imagesc(X0,Y0,I);
title('Input beam')
f1.Position(1:4) = [320 100 840 630];
xticks(-L0:L0/5:L0);yticks(-L0:L0/5:L0)
xlim([-L0 L0]);ylim([-L0 L0])
fontsize(1,25,"points")
axis tight;axis equal; axis off;
colormap(color)
hold on

%% Tightly focused image

%%% x-y　coordinate
N=501;
L=1.5;
%L=35; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi,r] = cart2pol(x,y);

%%% Parameter
NA=0.95;
f=1*10^3;
w0=0.61*lam/NA
zr=pi*w0^2/lam
z=zr.*0;
alpha=asin(NA); %Max angle

%%% Debye integral
Ex=0;Ey=0;Ez=0;
for t=0:60
    w=alpha/60;
    thita=w*t;
    x_int=f.*cos(phi).*tan(thita);
    y_int=f.*sin(phi).*tan(thita);
    ex=HGmode(n1,m1,x_int,y_int,w1,0,10^10,k).*ex_p;
    ey=HGmode(n2,m2,x_int,y_int,w1,0,10^10,k).*ey_p;
    J0=besselj(0,k.*r.*sin(thita)); J1=besselj(1,k.*r.*sin(thita)); J2=besselj(2,k.*r.*sin(thita));
    ax=sin(thita).*sqrt(cos(thita)).*exp(1j.*k.*z.*cos(thita)).*(ex.*J0.*(cos(thita)+1)+J2.*(1-cos(thita)).*(ex.*cos(2.*phi)+ey.*sin(2.*phi)));
    ay=sin(thita).*sqrt(cos(thita)).*exp(1j.*k.*z.*cos(thita)).*(ey.*J0.*(cos(thita)+1)+J2.*(1-cos(thita)).*(ex.*sin(2.*phi)-ey.*cos(2.*phi)));
    az=sin(thita).^2.*sqrt(cos(thita)).*exp(1j.*k.*z.*cos(thita)).*(1j.*J1.*(ex.*cos(phi)+ey.*sin(phi)));
    Ex=Ex+(ax).*w;
    Ey=Ey+(ay).*w;
    Ez=Ez+(az).*w;
end
Ex=1j.*Ex;
Ey=1j.*Ey;
Ez=2.*1j.*Ez;

%%% intensity plot
i_x=Ex.*conj(Ex);
i_y=Ey.*conj(Ey);
i_z=Ez.*conj(Ez);
i_t=i_x+i_y+i_z;

I_x=i_x./max(max(abs(i_t))); 
I_y=i_y./max(max(abs(i_t)));
I_z=i_z./max(max(abs(i_t)));
I_t=i_t./max(max(abs(i_t)));

%%% phase plot
i_xa=angle(i_x);
i_ya=angle(i_y);
i_za=angle(i_z);


%% Polarization strip
EX=Ex;
EY=Ey;
EZ=Ez;
EE=dot(EX,EX)+dot(EY,EY)+dot(EZ,EZ);

ax=real(EX.*sqrt(EE))./abs(sqrt(EE));
ax_n=ax./max(max(ax));
ay=real(EY.*sqrt(EE))./abs(sqrt(EE));
ay_n=ay./max(max(ay));
az=real(EZ.*sqrt(EE))./abs(sqrt(EE));
az_n=az./max(max(az));


cnt1=0;
cnt2=0;
g1=[];
g2=[];
r_int=(N-1)/2*(0.2/0.75);
for i=1:N
    for j=1:N
        I=i-(N+1)/2;
        J=j-(N+1)/2;
        if r_int*0.96<sqrt(I^2+J^2)
            if sqrt(I^2+J^2)<r_int*1.04
                cnt1=cnt1+1;
                g1(cnt1)=i+j*N;
                aX(cnt1)=ax(i,j);
                aY(cnt1)=ay(i,j);
                aZ(cnt1)=az(i,j);
                x_1(cnt1)=x(i,j);
                y_1(cnt1)=y(i,j);
                z_1(cnt1)=0;
            end
        end
    end
end


%% Figure show
f2=figure(2);
imagesc(X,Y,I_x);
hold on
title('ex')
%scatter(x_1,y_1,100,'white','.');
%imagesc(X,Y,i_xa);
f2.Position(1:4) = [320 100 840 630];
fontsize(2,25,"points")
xlim([-1000 1000]);ylim([-1300 1300])
xticks(-L:L/5:L);yticks(-L:L/5:L)
colormap(color)
colorbar
axis tight;axis equal; axis off;
hold off


f3=figure(3);
imagesc(X,Y,I_y);
hold on
title('ey')
%scatter(x_1,y_1,100,'white','.');
%imagesc(X,Y,i_ya);
f3.Position(1:4) = [320 100 840 630];
fontsize(3,25,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal; axis off;
hold off

f4=figure(4);
imagesc(X,Y,I_z);
hold on
title('ez')
%scatter(x_1,y_1,100,'white','.');
%imagesc(X,Y,i_za);
f4.Position(1:4) = [320 100 840 630];
fontsize(4,25,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal; axis off;
hold off

f5=figure(5);
imagesc(X,Y,I_t);
title('e')
%imagesc(X,Y,i_xa+i_ya+i_za);
f5.Position(1:4) = [320 100 840 630];
fontsize(5,25,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal; axis off;
hold on
%q=scatter(x_1,y_1,100,'white','.');
%q.MarkerFaceColor = [1 1 1];


f6=figure(6);
Linew=1;
scale=0.2;
q1=quiver3(x_1,y_1,z_1,aX.*scale,aY.*scale,aZ.*scale./2,1);
q1.LineWidth=Linew;
q1.Color=[1 0 0];
hold on
q2=quiver3(x_1,y_1,z_1,-aX.*scale,-aY.*scale,-aZ.*scale./2,1);
q2.LineWidth=Linew;
q2.Color=[0 0 1];
fontsize(6,15,"points")
title('polarization strip')

saveas(f1,'fig/input.png');saveas(f2,'fig/ex.png');saveas(f3,'fig/ey.png');saveas(f4,'fig/ez.png');saveas(f5,'fig/e.png');

z_plane=ones(1,numel(x_1)).*-0.05;
q3=quiver3(x_1,y_1,z_plane,aX.*scale,aY.*scale,aZ.*0);
q3.LineWidth=Linew;
q3.Color=[1 0.5 0.5];
q4=quiver3(x_1,y_1,z_plane,-aX.*scale,-aY.*scale,aZ.*0,'b');
q4.LineWidth=Linew;
q4.Color=[0.5 0.5 1];
xlim([-0.4 0.4])
ylim([-0.4 0.4])
zlim([-0.05 0.05])



%% Function
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
y=LG1./max(max(abs(LG1)));
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


function y=HGmode(n,m,x,y,w,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    pol=HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y);
    u=(1/W).*exp((-1).*(x.^2+y.^2)./W.^2).*exp(-(1j*k.*(x.^2+y.^2)./R)).*pol;%.*exp(gouy);
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
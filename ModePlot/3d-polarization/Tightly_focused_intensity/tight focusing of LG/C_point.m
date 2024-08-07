clear all
close all
clc

color='jet';
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%% Topological charge
l1=1; %North pole
l2=-1; %South pole1
%%% Radial index
p1=0;
p2=0;

%%% Polarizer angle 
ang=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of LG mode
w1=1000; %Beam waist
w2=w1;
lam=1.064; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2 %Rayleigh length

%z=input('position: '); %Beam position

%% Input beam image
%%% x-y　coordinate
N0=1001;
L0=1500;
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

ix_i=(LG_i1+LG_i2); 
iy_i=1j.*(LG_i1-LG_i2);
Ix_i=real(ix_i.*conj(ix_i));
Iy_i=real(iy_i.*conj(iy_i));
I=Ix_i./max(max(Ix_i))+Iy_i./max(max(Iy_i));%+Iz./max(max(Iz)); %Total input field intensity


f1=figure(1);
imagesc(X0,Y0,Iy_i);
title('Input beam')
f1.Position(1:4) = [320 100 840 630];
xticks(-L0:L0/5:L0);yticks(-L0:L0/5:L0)
xlim([-L0 L0]);ylim([-L0 L0])
fontsize(1,15,"points")
axis tight;axis equal;
colormap(color)
colorbar
hold on

%% Tightly focused image

%%% x-y　coordinate
N=151;
L=2;
%L=35; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi,r] = cart2pol(x,y);

%%% Parameter
NA=0.9;
f=1.7*10^3;
w0=10^3;
zr=pi*w0^2/lam
z=zr.*0;
alpha=asin(NA);
thitaP=phi;


Ex=1j;
Ey=1;

eX_x1=Ex.*i^l1.*(I0(l1,l1,alpha,k,z,phi,r,f)+I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX_x2=Ex.*i^l2.*(I0(l2,l2,alpha,k,z,phi,r,f)+I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX_y1=Ey.*i^(l1+1).*(-I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX_y2=Ey.*i^(l2+1).*(-I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX=(eX_x1+eX_x2)+(eX_y1-eX_y2);

eY_x1=Ex.*1j^(l1+1).*(-I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY_x2=Ex.*1j^(l2+1).*(-I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY_y1=Ey.*1j^l1.*(I0(l1,l1,alpha,k,z,phi,r,f)-I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)-I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY_y2=Ey.*1j^l2.*(I0(l2,l2,alpha,k,z,phi,r,f)-I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)-I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY=(eY_x1+eY_x2)+(eY_y1-eY_y2);

eZ_x1=Ex.*1j.^(l1-1).*(I2(l1,l1+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l1,l1-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ_x2=Ex.*1j.^(l2-1).*(I2(l2,l2+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l2,l2-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ_y1=Ey.*1j.^l1.*(-I2(l1,l1+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l1,l1-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ_y2=Ey.*1j.^l2.*(-I2(l2,l2+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l2,l2-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ=(eZ_x1+eZ_x2)+(eZ_y1-eZ_y2);

% figure(2)
% imagesc(X,-Y,angle(eZ))
% colormap('hsv')
% colorbar()
iX=eX.*conj(eX);
iY=eY.*conj(eY);
iZ=eZ.*conj(eZ);
it=iX+iY+iZ;
%iZ=real(eZ).^2;
IX=iX./max(max(iX));
IY=iY./max(max(iY));
IZ=iZ./max(max(iZ));
%%% phase plot
i_za=angle(eX);
i_xa=angle(eY);
i_ya=angle(eZ);


%% Polarization strip

EX=eX;
EY=eY;
EZ=eZ;
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
imagesc(X,-Y,iX);
hold on
title('ex')
scatter(x_1,y_1,100,'white','.');
%imagesc(X,Y,i_xa);
f2.Position(1:4) = [320 100 840 630];
fontsize(2,15,"points")
xlim([-1000 1000]);ylim([-1300 1300])
xticks(-L:L/5:L);yticks(-L:L/5:L)
colormap(color)
colorbar
axis tight;axis equal;
hold off


f3=figure(3);
imagesc(X,-Y,iY);
hold on
title('ey')
scatter(x_1,y_1,100,'white','.');
%imagesc(X,Y,i_ya);
f3.Position(1:4) = [320 100 840 630];
fontsize(3,15,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal;
hold off

f4=figure(4);
imagesc(X,-Y,iZ);
hold on
title('ez')
scatter(x_1,y_1,100,'white','.');
%imagesc(X,Y,i_za);
f4.Position(1:4) = [320 100 840 630];
fontsize(4,15,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal;
hold off

f5=figure(5);
imagesc(X,-Y,it);
title('e')
%imagesc(X,Y,i_xa+i_ya+i_za);
f5.Position(1:4) = [320 100 840 630];
fontsize(5,15,"points")
xticks(-L:L/5:L);yticks(-L:L/5:L)
xlim([-1300 1300]);ylim([-1300 1300])
colormap(color)
colorbar
shading interp; lighting phong; view(2); axis tight;axis equal;
hold on
q=scatter(x_1,y_1,100,'white','.');
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

z_adj=1;
z_plane=ones(1,numel(x_1)).*-0.05;
q3=quiver3(x_1,y_1,z_plane,aX./z_adj,aY./z_adj,aZ.*0);
q3.LineWidth=Linew;
q3.Color=[1 0.5 0.5];
q4=quiver3(x_1,y_1,z_plane,-aX./z_adj,-aY./z_adj,aZ.*0,'b');
q4.LineWidth=Linew;
q4.Color=[0.5 0.5 1];
xrange=1;
yrange=1;
zrange=1;
xlim([-xrange xrange])
ylim([-yrange yrange])
zlim([-zrange zrange])


%%%%%%%%%%%%%%Function%%%%%%%%%%%%%%%%%%%%
function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2./z; %Beam curvature
W=w.*(1+(z./zr).^2).^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
%Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*(r.^2./(2*R))+k.*z);
Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*z);
LG1=C./W.*exp((-1).*r.^2./W.^2).*(2^0.5*r./W).^abs(l).*laguerre(p, l, 2*r.^2./W.^2).*exp(Gouy);
LG=abs(LG1.^2);
%y=LG./max(LG(:));
y=LG1;

end

%% Function


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

function y=I0(l1,l2,alpha,k,z,phi,r,f)
    I=0;
    w0=1000;lam=1.064;
    for t=0:100
        w=alpha/100;
        thita=w*t;
        r_int=f.*tan(thita);
        A=LGmode(0,l1,r_int,phi,z,w0,lam);
        I=I+sqrt(cos(thita)).*sin(thita).*(1+cos(thita)).*besselj(l2,k.*r.*sin(thita)).*exp(1j.*k.*z.*cos(thita)).*A;
    end
    y=-1j*k.*f.*I;
end

function y=I1(l1,l2,alpha,k,z,phi,r,f)
    I=0;
    w0=1000;lam=1.064;
    for t=0:100
        w=alpha/100;
        thita=w*t;
        r_int=f.*tan(thita);
        A=LGmode(0,l1,r_int,phi,z,w0,lam);
        I=I+sqrt(cos(thita)).*sin(thita).*(1-cos(thita)).*besselj(l2,k.*r.*sin(thita)).*exp(1j.*k.*z.*cos(thita)).*A;
    end
    y=-1j*k.*f.*I/2;
end

function y=I2(l1,l2,alpha,k,z,phi,r,f)
    I=0;
    w0=1000;lam=1.064;
    for t=0:100
        w=alpha/100;
        thita=w*t;
        r_int=f.*tan(thita);
        A=LGmode(0,l1,r_int,phi,z,w0,lam);
        I=I+sqrt(cos(thita)).*sin(thita).^2.*besselj(l2,k.*r.*sin(thita)).*exp(1j.*k.*z.*cos(thita)).*A;
    end
    y=-1j*k.*f.*I;
end

function y=HGmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-(r.^2./W^2)-(1j*k.*(x.^2+y.^2)./R));%.*exp(gouy);
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
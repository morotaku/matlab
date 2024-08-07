clear all
close all
clc
%%% Topological charge
l1=0;
l2=-1;
%%% Radial index
p=0;
scale=7;
phi=0;
thita=pi/2;
%%% Polarizer angle 
angle=input('polarizer angle: '); %thita=-1 → no polarizerG
%%% Parameter
w=5; %Beam waist
lam=0.64; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size


%%% x-y　coordinate
N=1000;
L=10; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
lg=LGmode(p,l1,r,phi1,z,w,lam);
LG=lg.*conj(lg);
LG_n=LG./max(max(LG));
% imagesc(LG_n);
% colormap('gray')
% shading interp; lighting phong; view(2); axis equal; axis tight; axis off;
% hold on



%%%%%%%%%%%%%%quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=9;
X2=linspace(-L.*0.8,L.*0.8,N2);
Y2=linspace(-L.*0.8,L.*0.8,N2);
[x2,y2]=meshgrid(X2,Y2);
%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);
%%% LG mode
lg1=LGmode(p,l1,r2,phi2,z,w,lam);
lg2=LGmode(p,l2,r2,phi2,z,w,lam);
f=1;
imagesc([-L L],[-L L],LG_n);
colormap('gray')
hold on
for t = 1:2

    %%%North pole of Poincare sphere
    Er=sin(thita/2).*exp(1j*phi/2).*exp(1j*t*pi).*lg1;  %ℓ=1
    %%%South pole
    El=cos(thita/2).*exp(-1j*phi/2).*exp(1j*t*pi).*lg2; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1=real(Er);
    ey1=real(-1j.*Er);
    %%%|RHC> ℓ=-1
    ex2=real(El);
    ey2=real(1j.*El);
    
    ex=ex1+ex2;
    ey=ey1+ey2;
    [Ex,Ey]=polarizer(angle,ex,ey);
    %Ex=Ex.*conj(Ex);Ey=Ey.*conj(Ey);
    q=quiver(x2,y2,Ex.*scale,Ey.*scale,'off');
    q.LineWidth=3;
    q.Color='red';
    axis equal; axis tight; axis off;
    hold on
    pause(0.1)
%     myMovie(f) = getframe(gcf);
%     f=f+1;
end

% v = VideoWriter('myMovie.mp4', 'MPEG-4');
% % ゆっくり表示したかったのでフレームレートを初期値より下げています
% v.FrameRate = 10;
% open(v);
% writeVideo(v, myMovie)
% close(v);

%imagesc(i_a);
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


clear all
close all
clc
%%% Topological charge
l1=1; %North pole
l2=-1; %South pole
%%% Radial index
p1=0;
p2=0;
%%% Poincare sphere angle
phi=pi; %latitude
thita=pi/2; %longitude

%%% x-y　coordinate
L=5;
N2=11;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
Z2=-12:4:8;
[x2,y2,z]=meshgrid(X2,Y2,Z2);


%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%



%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

t=0;
%%% Parameters
w=1; %Beam waist
lam=0.64; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2./z; %Beam curvature
R(:,:,1)=R(:,:,1)*0.1^10;
W=w*(1+(z./zr).^2).^0.5; %Beam size

wave_angle=abs(asin(r2./R));
wave_angle(:,:,1)=pi/2;
%%% Polarizer angle
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer


%%% LG mode
lg1_q=LGmode(p1,l1,r2,phi2,z,w,lam);
lg2_q=LGmode(p1,l2,r2,phi2,z,w,lam);
LG1_q=LGmode(p1,l1,r2,phi2,z,w,lam)./max(max(lg1_q));
LG2_q=LGmode(p1,l2,r2,phi2,z,w,lam)./max(max(lg1_q));
%%%North pole of Poincare sphere
E1_q=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2)).*exp(1j*t*pi/20).*LG1_q;
%%%South pole
E2_q=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2)).*exp(1j*t*pi/20).*LG2_q; %ℓ=-1
%%%|LHC> ℓ=1
ex1_q=real(E1_q);
ey1_q=real(-1j.*E1_q);
%%%|RHC> ℓ=-1
ex2_q=real(E2_q);
ey2_q=real(1j.*E2_q);

ex_q=ex1_q+ex2_q;
ey_q=ey1_q+ey2_q;
%[Ex_q,Ey_q]=polarizer(angle,ex_q,ey_q);
Ex_q=ex_q./max(max(ex_q(:)));
Ey_q=ey_q./max(max(ey_q(:)));
for i=1:numel(Ex_q)
    if Ex_q(i)>0.3
    elseif Ex_q(i)<-0.3
    else
        Ex_q(i)=0;
    end
end
for i=1:numel(Ey_q)
    if Ey_q(i)>0.3
    elseif Ey_q(i)<-0.3
    else
        Ey_q(i)=0;
    end
end
Ez_q=zeros(size(z));%LG1_q.*abs(cos(wave_angle));


    
    
zero=zeros(N2,N2,numel(Z2));
q=quiver3(x2,z,y2,Ex_q,Ez_q,Ey_q,0);
q.LineWidth=1.75;
q.Color='red';
set(gca,'TickLength',[0 0]);
axis equal;%axis off



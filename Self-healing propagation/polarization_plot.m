clear all
close all
clc

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude強度比
ang=-1;
L=300;
z=10.*10^4;
%%%Figure parameter
width=10;
scale=L/16;
%%%%%%%%%%%%%%Intensity%%%%%%%%%%%%%%%%%%

N1=300;
%L=35; %Display lange
X=linspace(-L,L,N1);
Y=linspace(-L,L,N1);


beam1_x=self_healing(z,0,N1,L);
beam1_y=self_healing(z,0,N1,L).*exp(1j.*pi/2);
beam2_x=self_healing(z,-1,N1,L);
beam2_y=self_healing(z,-1,N1,L).*exp(-1j.*pi/2);
beam_x=beam1_x+beam2_x;
beam_y=beam1_y+beam2_y;
Beam=abs(beam_x.*conj(beam_x))+abs(beam_y.*conj(beam_y));
Beam_n=Beam./max(max(Beam));


imagesc(X,Y,Beam);
colormap('gray');
%colormap('jet');
axis equal; 
shading interp; lighting phong; view(2); axis tight;
hold on

%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=19;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);
%%% LG mode
LG_n1=self_healing(z,0,N2,L);

%%% HGB
LG_n2=self_healing(z,-1,N2,L);



cnt1=0;
cnt2=0;
g1=[];
g2=[];

for i=1:N2^2
    if cos(thita/2)*abs(LG_n1(i))>sin(thita/2)*abs(LG_n2(i))
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
    [Ex_1,Ey_1]=polarizer(ang,ex_1,ey_1);
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
    [Ex_2,Ey_2]=polarizer(ang,ex_2,ey_2);
    Ex2=Ex_2.*scale+x_2;
    Ey2=Ey_2.*scale+y_2;
    axis equal;axis tight; axis off;
    
    scatter(Ex2,Ey2,width,'.','blue')
    pause(0.005)
end

% %%%Right circular
% for t = 0:120
%     %%%North pole of Poincare sphere
%     Er_1=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/60).*LG_r1;  %ℓ=1
%     %%%South pole
%     El_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*LG_l1; %ℓ=-1
%     %%%|LHC> ℓ=1
%     exr_1=Er_1;
%     eyr_1=1j.*Er_1;
%     %%%|RHC> ℓ=-1
%     exl_1=El_1;
%     eyl_1=-1j.*El_1;
%     
%     ex_1=(exr_1+exl_1).*conj(exr_1+exl_1);
%     ey_1=(eyr_1+eyl_1).*conj(eyr_1+eyl_1);
%     [Ex_1,Ey_1]=polarizer(ang,ex_1,ey_1);
%     Ex1=Ex_1.*scale+x_1;
%     Ey1=Ey_1.*scale+y_1;
%     axis equal; 
%     
%     scatter(Ex1,Ey1,width,'.','red')
%     pause(0.005)
% end
% hold on
% fig=gcf;
% fig.Color = 'none';
% fig.InvertHardcopy = 'off';
% 
% %%%Left circular
% for t = 0:120
%     %%%North pole of Poincare sphere
%     Er_2=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/60).*LG_r2;%.*exp(1j.*(l1.*phi2(g2)))
%     %%%South pole
%     El_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*LG_l2; %.*exp(1j.*(l2.*phi2(g2)))
%     %%%|LHC> ℓ=1
%     ex1_2=Er_2;
%     ey1_2=1j.*Er_2;
%     %%%|RHC> ℓ=-1
%     ex2_2=El_2;
%     ey2_2=-1j.*El_2;
%     
%     ex_2=(ex1_2+ex2_2).*conj(ex1_2+ex2_2);
%     ey_2=(ey1_2+ey2_2).*conj(ey1_2+ey2_2);
%     [Ex_2,Ey_2]=polarizer(ang,ex_2,ey_2);
%     Ex2=Ex_2.*scale+x_2;
%     Ey2=Ey_2.*scale+y_2;
%     axis equal;axis tight; axis off;
%     
%     scatter(Ex2,Ey2,width,'.','blue')
%     pause(0.005)
% end


function y=self_healing(z,l,N,L)
b=0;
w0=45; % obstacle
lam=1; % wavelength
angle_inc=0.5*pi/180; %    radians
n=1.5; %  refractive index
k=2*pi/lam;
kt=k.*(n-1).*sqrt(1-cos(angle_inc).^2);

%x-y　coordinate
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[thita,r] = cart2pol(x,y);
%thita=pi/6;
beta=0;

%phi=atan(2.*b.*z.*sin(beta)-1j.*k.*r.*w0^2.*sin(thita)./(2.*b.*z.*cos(beta)-1j.*k.*r.*w0^2.*cos(thita)));
phi=atan2(sin(thita),cos(thita));
A=sqrt(k^2.*r.^2./z.^2-4.*b.^2/w0.^4+4.*1j.*b.*k.*r.*cos(thita-beta)./(w0^2.*z));


u2_amp=k*w0^2./(k*w0^2+2.*1j.*z).*besseli(l,1j.*kt.*A.*w0^2.*z./(k.*w0^2+2.*1j.*z)).*exp(1j.*l.*(pi/2-phi));
u2_phz=exp(1j.*k.*z.*(1+r.^2/(2.*z.^2)-(kt.^2+A.^2).*w0^2./(2.*k^2*w0^2+4.*1j.*k.*z))-b.^2./w0.^2);
u2=u2_amp.*u2_phz;
u1=exp(1j.*k.*(z-kt^2/(2*k^2))+1j.*l.*(pi-thita)).*besselj(l,kt.*r);

u=u1-u2;
% I=abs(u.*conj(u));
% I_n=I./max(max(I));
y=u./max(max(abs(u)));
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


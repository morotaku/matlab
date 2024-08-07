clear all
close all
clc
%% Parameters
%Pump_beam_index
l_p=0;
p_p=0;

%Oscillating_beam_index
l_l=0;
p_l=0;

%Pump_beam_Parameters
w_p=50; %Beam size
lam_p=0.442; %Wave length
z_p=0; %Beam position 

%Oscillating_beam_Parameters
w_l=10; %Beam size
lam_l=0.640; %Wave length
z_l=0; %Beam position 

%x-y　coordinate
N=1000;
L=500; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);

%Pump_Beam_spatial_profile
pump=LGmode(p_p,l_p,r,phi,z_p,w_p,lam_p);
Pump=abs(pump.*conj(pump));
Pump_n=Pump./max(max(Pump));

%Lasing_beam_spatial_profile
cav=LGmode(p_l,l_l,r,phi,z_l,w_l,lam_l);
Cav=abs(cav.*conj(cav));
Cav_n=Cav./max(max(Cav));

%% Spatial_profile(In case of Gaussian beam)
figure(1)
subplot(1,2,1);imagesc(Pump_n);
axis xy; axis equal; axis tight; %axis off;
colormap('hot')
subplot(1,2,2);imagesc(Cav_n);
axis xy; axis equal; axis tight; %axis off;
colormap('hot')
title('Spatial profile (l=0)')

%% Overlap_efficiency(l:Variable, w_l:Constant)
l_lv=linspace(1,30,30);
for i=1:30
    %Lasing_beam_spatial_profile
    cav=LGmode(p_l,l_lv(i),r,phi,z_l,w_l,lam_l);
    Cav=abs(cav.*conj(cav));
    Cav_n=Cav./max(max(Cav));
    E_eff(i)=dot(Pump_n(:),Cav_n(:)).^2./(dot(Pump_n(:),Pump_n(:)).*dot(Cav_n(:),Cav_n(:)));
end
figure(2)
plot(l_lv,E_eff)
xlabel('Topological charge');ylabel('Overlap Efficiency')



%% Overlap_efficiency(l:Constant, w_l:Variable)
w_l=linspace(10,100,30);
for i=1:30
    %Lasing_beam_spatial_profile
    cav=LGmode(p_l,l_l,r,phi,z_l,w_l(i),lam_l);
    Cav=abs(cav.*conj(cav));
    Cav_n=Cav./max(max(Cav));
    E_eff(i)=dot(Pump_n(:),Cav_n(:)).^2./(dot(Pump_n(:),Pump_n(:)).*dot(Cav_n(:),Cav_n(:)));
end
figure(3)
plot(w_l,E_eff)
xlabel('Lasing beam size');ylabel('Overlap Efficiency')







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

function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2./z; %Beam curvature
W=w.*(1+(z./zr).^2).^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*(r.^2./(2*R))+k.*z);
LG=C./W.*exp((-1).*r.^2./W.^2).*(2^0.5*r./W).^abs(l).*laguerre(p, l, 2*r.^2./W.^2).*exp(Gouy);
y=LG;
end

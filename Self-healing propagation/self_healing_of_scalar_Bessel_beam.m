close all
clear all
clc

figure(1)
colormap('jet')
z_0=transpose(self_healing(0.000001));
imagesc(z_0);
shading interp; lighting phong;  axis equal; axis tight; axis xy

N=200;
Z=linspace(0.*10^4,15.*10^4,N);
for i=1:N
    z=Z(i);
    I=self_healing(z);
    propagation(i,1:N)=I(N/2,:);
end
figure(2)
colormap('jet')
propagation_t=transpose(propagation);
imagesc(propagation_t);
shading interp; lighting phong;  axis equal; axis tight; axis xy

function y=self_healing(z)
l=2;
b=200;
w0=100; % obstacle
lam=1; % wavelength
angle_inc=0.5*pi/180; %    radians
n=1.5; %  refractive index
k=2*pi/lam;
kt=k.*(n-1).*sqrt(1-cos(angle_inc).^2);

%x-y　coordinate
N=200;
L=400; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[thita,r] = cart2pol(x,y);
%thita=pi/6;
beta=0;

phi=atan(2.*b.*z.*sin(beta)-1j.*k.*r.*w0^2.*sin(thita)./(2.*b.*z.*cos(beta)-1j.*k.*r.*w0^2.*cos(thita)));
A=sqrt(k^2.*r.^2./z.^2-4.*b.^2/w0.^4+4.*1j.*b.*k.*r.*cos(thita-beta)./(w0^2.*z));


u2_amp=k*w0^2./(k*w0^2+2.*1j.*z).*besseli(l,1j.*kt.*A.*w0^2.*z./(k.*w0^2+2.*1j.*z)).*exp(1j.*l.*(pi/2-phi));
u2_phz=exp(1j.*k.*z.*(1+r.^2/(2.*z.^2)-(kt.^2+A.^2).*w0^2./(2.*k^2*w0^2+4.*1j.*k.*z))-b.^2./w0.^2);
u2=u2_amp.*u2_phz;
u1=exp(1j.*k.*(z-kt^2/(2*k^2))+1j.*l.*(pi-thita)).*besselj(l,kt.*r);

u=u1-u2;
I=abs(u.*conj(u));
I_n=I./max(max(I));
y=I_n;
end


% phi=atan(2.*b.*z.*sin(beta)-1j.*k.*r.*w0^2.*sin(thita)./(2.*b.*z.*cos(beta)-1j.*k.*r.*w0^2.*cos(thita)));
% A=sqrt(k^2.*r.^2./z.^2-4.*b.^2/w0.^4+4.*1j.*b.*k.*r.*cos(thita-beta)./(w0^2.*z));
% 
% 
% u2_amp=k*w0^2/(k*w0^2+2.*1j.*z).*besseli(l,1j.*kt.*A.*w0^2.*z/(k.*w0^2+2.*1j.*z)).*exp(1j.*l.*(pi/2-phi));
% u2_phz=exp(1j.*k.*z.*(1+r.^2/(2.*z^2)-(kt.^2+A.^2).*w0^2./(2.*k^2*w0^2+4.*1j.*k.*z))-b.^2./w0.^2);
% u2=u2_amp.*u2_phz;
% u1=exp(1j.*k.*(z-kt^2/(2*k^2))+1j.*l.*(pi-thita)).*besselj(l,kt.*r);
% 
% u=u2-u1;
% I=abs(u).^2;
% I_n=I./max(max(I));


%modified bessel function besseli(n,x)
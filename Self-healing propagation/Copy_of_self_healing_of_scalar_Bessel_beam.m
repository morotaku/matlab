close all
clear all
clc

figure(1)
colormap('jet')
beam1_x=self_healing(0.00000001,0);
beam1_y=self_healing(0.00000001,0).*exp(1j.*pi/2);
beam2_x=self_healing(0.00000001,-1);
beam2_y=self_healing(0.00000001,-1).*exp(-1j.*pi/2);
beam_x=beam1_x+beam2_x;
beam_y=beam1_y+beam2_y;
Beam=abs(beam_x.*conj(beam_x))+abs(beam_y.*conj(beam_y));
Beam_n=Beam./max(max(Beam));


imagesc(Beam);
shading interp; lighting phong;  axis equal; axis tight; axis xy

N=250;
L=300;
X=linspace(-L,L,N);
Z=linspace(0.*10^4,15.*10^4,N);
for i=1:N
    z=Z(i);
    beam1_x=self_healing(z,0);
    beam1_y=self_healing(z,0).*exp(1j.*pi/2);
    beam2_x=self_healing(z,-1);
    beam2_y=self_healing(z,-1).*exp(-1j.*pi/2);
    beam_x=beam1_x+beam2_x;
    beam_y=beam1_y+beam2_y;
    Beam=abs(beam_x.*conj(beam_x))+abs(beam_y.*conj(beam_y));
    Beam_n=Beam./max(max(Beam));
    propagation(i,1:N)=Beam_n(N/2,:);
%     if i==N
%         figure(2)
%         imagesc(Beam_n);
%         colormap('jet')
%         shading interp; lighting phong;  axis equal; axis tight; axis xy
%     end
end
figure(2)
colormap('jet')
propagation_t=transpose(propagation);
imagesc(Z,X,propagation_t);
axis tight; axis xy

function y=self_healing(z,l)
b=200;
w0=100; % obstacle
lam=1; % wavelength
angle_inc=0.5*pi/180; %    radians
n=1.5; %  refractive index
k=2*pi/lam;
kt=k.*(n-1).*sqrt(1-cos(angle_inc).^2)

%x-y　coordinate
N=250;
L=300; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[thita,r] = cart2pol(x,y);
%thita=pi/6;
beta=0;

phi=atan(2.*b.*z.*sin(beta)-1j.*k.*r.*w0^2.*sin(thita)./(2.*b.*z.*cos(beta)-1j.*k.*r.*w0^2.*cos(thita)));
%phi=atan2(sin(thita),cos(thita));

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
clear all
close all
clc
Red=[1 1 1];
for i=0:50
    Red=[Red;[1 1-0.02*i 1-0.02*i]];
end
for i=0:10
    Red=[Red;1 0 0];
end

%Radial index
l=0;

%Azimuth index
p=0;

%Parameters(um)
w=5; %Beam waist
lam=0.64; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
z=input('position: '); %Beam position 
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size

%x-y　coordinate
N=1000;
L=10; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);

lg=LGMODE(p,l,r,phi,z,w,lam);
LG=lg.*conj(lg);
LG_n=LG./max(max(LG));
LG_a=angle(lg);

figure(1)
imagesc(LG_n);
shading interp; lighting phong; view(2); axis tight;axis equal; axis off;
colormap('hot')
colorbar('Ticks',[0.01,1],...
         'TickLabels',{'MIN','MAX'});
fontsize(1,30,"points")

figure(2)
imagesc(LG_a);
shading interp; lighting phong; view(2); axis tight;axis equal; axis off;
colormap('gray');
colorbar('Ticks',[-3,3],...
         'TickLabels',{'0','2π'})
fontsize(2,40,"points")

function y=LGMODE(p,l,r,phi,z,w,lam)
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



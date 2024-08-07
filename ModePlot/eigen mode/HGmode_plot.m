clear all
close all
clc
%Mode index
n=0;
m=2;

%Parameters
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
%Polar coodinate
[phi,r] = cart2pol(x,y);

%hg=HGmode(n,m,x,y,w,r,z,R,k);
hg=2.*HGmode(1,0,x,y,w,r,z,R,k)+HGmode(0,2,x,y,w,r,z,R,k);
HG=hg.*conj(hg);
HG_n=HG./max(max(HG));
HG_a=angle(hg);

figure(1)
imagesc(HG_n)
colormap('jet')
shading interp; lighting phong; view(2); axis equal; axis tight; axis off;

figure(2)
imagesc(HG_a)
colormap('jet')
shading interp; lighting phong; view(2); axis equal; axis tight; axis off;
colorbar()

function y=HGmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*exp(-(r.^2./W^2)-(1j*k.*r.^2/(2*R))).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));
end

function u=HGmode45(n,x,y,r,w,z,zr,k,R)
    HG=0;
    for i=0:n
        HG=HG+nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
    end
    u=HG.*sqrt(2)^n;
end

function u=HGmode_45(n,x,y,r,w,z,zr,k,R)
    HG=0;
    for i=0:n
        if rem(i,2)==0
            HG=HG-nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        else
            HG=HG+nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        end
    end
    u=HG.*sqrt(2)^n;
end

function y=HermitePol(n,x)
Hn=0;
for m=0:floor(n/2)
    Hn=Hn+((-1)^m.*(2.*x).^(n-2*m))./(factorial(m)*factorial(n-2*m));
end
Hn=Hn.*factorial(n);
y=Hn;
end


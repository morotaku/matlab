%%% Simple computer program to calculate arbitrary tightly focused (propagating andevanescent) vector light fields
clear all; close all; clc;
n1=0;m1=2;
n2=1;m2=0;

L=3000; % main domain size
range=80000; %display range
M=max([n1 n2 m1 m2])+1;
crop=L/70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%initial parameters
z0=0; %interface position (focus at z=0)
z=linspace(-2,2,41); %axial domain
m0=1/2; %coordinate shift
xmax=10; %half the size of the cropped output in wavelengths
%(final size of images is 2xmax)
ref1=1; % refractive index for incident field (immersion oil)
ref2=1; % refractive index for reflected field (immersion oil)
ref3=1; % refractive index for transmitted field (after interface)
e1=ref2^2;
e2=ref3^2;
f=1800; % effective focal length of the lens in um (10ˆ−6 m)
NA=0.95; %effective numerical aperture
lam=0.64; %wavelength in um
k=2*pi/lam; %Wave number
l=0.64; %wavelength u
R=f*NA/ref2; %aperture radius in um;
w0=f*NA/ref2; % The beam waist is equal to the aperture radius.
N=w0.*1.*(sqrt(M)); % number of points for the aperture radius
dx_inf=R/N; %spatial dimension at aperture um/px or um/point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Spatial coordinate system before lens (x inf, y inf)
d_inf=linspace(-range/2,range/2,L)*dx_inf; %spatial axis shifted by m0.
[x_inf,y_inf]=meshgrid(d_inf,-d_inf); %mesh, y inverted in Matlab images
[theta_inf,rho_inf] = cart2pol(x_inf,y_inf); %auxiliary polar coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% angular
k0=2*pi/lam;
k1=ref2*k0;
k2=ref3*k0;
dk=k0*NA/N; %dk=(1/2)*(k1/(f))*dx0;
kx=x_inf*k1/f;
ky=y_inf*k1/f;
dxf=2*pi/(L*dk); %equivalent to dxf=N*lambda/(L*NA);
%conversion factor at output plane in um/px
kM1=k1*ones(L,L); %%%array with magnitude of k1
kz1=sqrt(kM1.^2-kx.^2-ky.^2);
kM2=k2*ones(L,L);
kz2=sqrt(kM2.^2-kx.^2-ky.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmax2=round(xmax*lam/dxf); %xmax in px at output plane
xmax2=110;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Phase
%%%%%%%%%%%%%%%%%%%%%%%correction for m0 shift
x=linspace(-range/2,range/2,L); %%%auxiliary coordinate for correction
y=x; [X,Y]=meshgrid(x,y);
[phi,r] = cart2pol(X,Y);

%correction phase
PhaseShift=m0*2*pi*X/L+m0*2*pi*Y/L+0.8; %correction phase
PhaseShift=PhaseShift.*0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input intensity
%%%%%%%%%%%%%%%%%%%%%input field initial field transverse

%e_incx=HGmode(1,0,X,Y,w0,0,10^10,k); %normalized amplitude E x0 A 0
e_incx=HGmode(n1,m1,X,Y,w0,0,R,k);
%e_incx=LGmode(0,1,r,phi,0,w0,lam);%+exp(1j.*0).*LGmode(0,0,r,phi,0,w0,lam);

%e_incy=-HGmode(0,1,X,Y,w0,0,10^10,k);%.*1j;
e_incy=HGmode(n2,m2,X,Y,w0,0,R,k);
%e_incy=-1j.*(LGmode(0,1,r,phi,0,w0,lam)-exp(1j.*0).*LGmode(0,0,r,phi,0,w0,lam));

phi_incx=angle(e_incx);
phi_incy=angle(e_incy);
E_incx=abs(e_incx); E_incy=abs(e_incy);
E_incx=E_incx./max(max(E_incx)); E_incy=E_incy./max(max(E_incy));

figure(1)
fontsize(1,25,"points");colormap('jet');
subplot(2,2,1);imagesc(E_incx);axis image; colorbar(); title('Input intensity(x)')
subplot(2,2,2);imagesc(E_incy);axis image; colorbar(); title('Input intensity(y)')
subplot(2,2,3);imagesc(phi_incx);%axis image;colorbar(); clim([-pi pi]); title('Input phase(x)') 
subplot(2,2,4);imagesc(phi_incy);%axis image;colorbar(); clim([-pi pi]); title('Input phase(y)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%round soft aperture (spatial coordinates)
Theta=0.5*(1+tanh((1.5/1)*(N-sqrt((x_inf/dx_inf).^2+(y_inf/dx_inf).^2))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%polarization
%cx=ones(L,L); cy=zeros(L,L);   % lin horizontal
%cx=zeros(L,L); cy=ones(L,L); % lin vertical
%cx=1*ones(L,L); cy=1*ones(L,L); % lin diag
cx=1i*ones(L,L); cy=1i*ones(L,L); % lin diag
%cx=ones(L,L); cy=1i*ones(L,L); %%right circ
%cx=1i*ones(L,L); cy=ones(L,L); %%left circ
%cx=(kx)./sqrt(kx.ˆ2+ky.ˆ2); cy=(ky)./sqrt(kx.ˆ2+ky.ˆ2); %%radial cx=cos(phi) cy=sin(phi)
%cx=−ky./sqrt(kx.ˆ2+ky.ˆ2); cy= kx./sqrt(kx.ˆ2+ky.ˆ2);%%azimuthal cx=cos(phi) cy=sin(phi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Factors for E inf
CF1=-(1./kz1).*sqrt(kz1./kM1).*(1./(kx.^2+ky.^2))*1i*f*exp(-1i*f*2*pi*ref2/lam); %Factors for E_inf_r and f
CF2=-(1./kz2).*sqrt(kz1./kM1).*(1./(kx.^2+ky.^2))*1i*f*exp(-1i*f*2*pi*ref2/lam).*(kz2./kz1); %Factors for E inf t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%Focused propagating before interface
%for E inf
% x component
MTx1=(ky.^2 +kx.^2.*kz1./kM1).*CF1; %lens transmision factor, incident x
MTx2=(-ky.*kx +(ky.*kx.*kz1)./kM1).*CF1; %lens transmission factor incident y
%y component
MTy1=(-kx.*ky+(kx.*ky.*kz1)./kM1).*CF1;%lens transmision factor, incident x
MTy2=(kx.^2 +ky.^2.*kz1./kM1).*CF1;%lens transmision factor, incident y
%z component
MTz1=(-(kx.^2+ky.^2).*kx./kM1).*CF1;%lens transmision factor, incident x
MTz2=(-(kx.^2+ky.^2).*ky./kM1).*CF1; %lens transmision factor, incident y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%start
z=0;
Hzf=exp(1i*kz1.*z).*Theta; %propagating

%%%%%%%%%%% focal
% x component E inf
Fieldxxf=(cx.*MTx1.*E_incx.*exp(1i*phi_incx)+cy.*MTx2.*E_incy.*exp(1i*phi_incy)).*Hzf;
%integration x component
Ex_inff=ifftshift(ifft2(fftshift((Fieldxxf))));%.*exp(1i*PhaseShift);
Exf=Ex_inff(round(L/2+1-crop):round(L/2+1+crop),round(L/2+1-crop):round(L/2+1+crop)); %cropped
%Exf=Ex_inff;

%y component E inf
Fieldxyf=(cx.*MTy1.*E_incx.*exp(1i*phi_incx)+cy.*MTy2.*E_incy.*exp(1i*phi_incy)).*Hzf;
%integration y component
Ey_inff=ifftshift(ifft2((fftshift(Fieldxyf)))).*exp(1i*PhaseShift);
Eyf=Ey_inff(round(L/2+1-crop):round(L/2+1+crop),round(L/2+1-crop):round(L/2+1+crop));
%Eyf=Ey_inff;

%z component E inf
Fieldxzf=(cx.*MTz1.*E_incx.*exp(1i*phi_incx)+cy.*MTz2.*E_incy.*exp(1i*phi_incy)).*Hzf;
%integration z component
Ez0f=ifftshift(ifft2((fftshift(Fieldxzf)))).*exp(1i*PhaseShift);
Ezf=Ez0f(round(L/2+1-crop):round(L/2+1+crop),round(L/2+1-crop):round(L/2+1+crop));



Ex=abs(Exf.*conj(Exf)); % total field at the material n 2 incident+reflected
Ey=abs(Eyf.*conj(Eyf));
Ez=abs(Ezf.*conj(Ezf));
Et=Ex+Ey+Ez;
et=max(max(Et));
Ex_n=Ex./max(max(et));
Ey_n=Ey./max(max(et));
Ez_n=Ez./max(max(et));
Et_n=Et./max(max(et));
%%%%%%%%%%%
figure(2)
subplot(1,3,1);imagesc(abs(Fieldxxf.*conj(Fieldxxf)));axis image; colorbar(); title('x');
subplot(1,3,2);imagesc(abs(Fieldxyf.*conj(Fieldxyf)));axis image; colorbar(); title('y')
subplot(1,3,3);imagesc(abs(Fieldxzf));axis image; colorbar(); title('z')
round(max(max(Ex_n)),2)

phase_x=angle(Exf);phase_y=angle(Eyf);phase_z=angle(Ezf);
%phase_x(phase_x>2.8)=-phase_x(phase_x>2.8); phase_y(phase_y>2.8)=-phase_y(phase_y>2.8); phase_z(phase_z>2.8)=-phase_z(phase_z>2.8);

figure(3)
colormap('hot')
subplot(1,4,1);imagesc(Et_n);axis image; colorbar;axis off; clim([0 1]); title('Total Intensity');
subplot(1,4,2);imagesc(Ex_n);axis image; colorbar;axis off; clim([0 max(max(Ex_n))]); title('x');
subplot(1,4,3);imagesc(Ey_n);axis image; colorbar;axis off; clim([0 max(max(Ey_n))]); title('y');
% subplot(1,4,2);imagesc(Ex_n);axis image; colorbar;axis off; clim([0 1]); title('x');
% subplot(1,4,3);imagesc(Ey_n);axis image; colorbar;axis off; clim([0 1]); title('y');
subplot(1,4,4);imagesc(Ez_n);axis image; colorbar;axis off; clim([0 max(max(Ez_n))]); title('z');
fontsize(3,22,"points");

figure(4)
colormap('hsv')
subplot(1,4,2);imagesc(phase_x);axis image; colorbar(); axis off; clim([-pi pi]);
subplot(1,4,3);imagesc(phase_y);axis image; colorbar(); axis off; clim([-pi pi]);
subplot(1,4,4);imagesc(phase_z);axis image; colorbar(); axis off; clim([-pi pi]);
fontsize(4,22,"points");




%% Function
function y=HGmode(n,m,x,y,w,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    %u=(1/W).*exp(-((x.^2+y.^2)./W^2)-(1j*k.*(x.^2+y.^2)./(2*R))).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y);%.*exp(gouy);
    u=(1/W).*exp(-((x.^2+y.^2)./W^2)).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y);%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));
end

function u=HGmode45(n,x,y,w,z,k,R)
    HG=0;
    for i=0:n
        HG=HG+nchoosek(n,i).*HGmode(i,n-i,x,y,w,z,R,k);
    end
    u=HG.*sqrt(2)^n;
end

function u=HGmode_45(n,x,y,w,z,k,R)
    HG=0;
    for i=0:n
        if rem(i,2)==0
            HG=HG-nchoosek(n,i).*HGmode(i,n-i,x,y,w,z,R,k);
        else
            HG=HG+nchoosek(n,i).*HGmode(i,n-i,x,y,w,z,R,k);
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

function y=rotate(img,ang)
    % 回転する角度を指定 (例: 45度)
    rotationAngle = ang;
    
    % 画像のサイズを取得
    [rows, cols, ~] = size(img);
    
    % 回転行列の作成
    theta = deg2rad(rotationAngle);
    rotationMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    % 中心座標の計算
    center = [(cols+1)/2, (rows+1)/2];
    
    % 回転行列を用いて画像を回転
    rotatedImage = zeros(size(img));
    for i = 1:rows
        for j = 1:cols
            % 回転前の座標を計算
            originalCoord = rotationMatrix * ([j; i] - center') + center';
            
            % 回転前の座標が画像の範囲内か確認してから値を代入
            if all(originalCoord > 0) && all(originalCoord <= [cols; rows])
                % バイリニア補間を行う場合は、適切な補間関数を実装してください
                
                rotatedImage(i, j, :) = img(ceil(originalCoord(2)), ceil(originalCoord(1)), :);
            end
        end
    end
    y=rotatedImage;
end
%%% Simple computer program to calculate arbitrary tightly focused (propagating andevanescent) vector light fields
clear all; close all; clc;
M=5;n1=1; n2=1;
N=1300; % number of points for the aperture radius
range=110000; %display range
L=4000;
ratio=1;
fig_size=200; %Inputビームの画像サイズ変更


%%%傾き補正
til1=1;
til2=1;
%%%切片補正
itc1=0;
itc2=0;
pointx_0=2200; %左端のローブ
pointx_n=2300;
pointy_0=2145; %左端のローブ
pointy_n=2370;


imp1=importdata('0.csv'); imp1_r=imresize(imp1,[fig_size fig_size]);
imp_x=abs(sqrt(imp1_r)); imp_x=imp_x./max(max(imp_x));


imp2=importdata('90.csv'); imp2_r=imresize(imp2,[fig_size fig_size]);
imp_y=abs(sqrt(imp2_r)); imp_y=imp_y./max(max(imp_y));

L2=size(imp_x);
L2=L2(1);

adj=length(L/2-L2/2:L/2+L2/2)-L2;
e_incx=zeros(L); e_incy=zeros(L);
e_incx(L/2-L2/2:L/2+L2/2-adj,L/2-L2/2:L/2+L2/2-adj)=sqrt(imp_x);
e_incy(L/2-L2/2:L/2+L2/2-adj,L/2-L2/2:L/2+L2/2-adj)=sqrt(imp_y);
% figure(1);imagesc(e_incx); axis equal; axis tight;
% figure(2);imagesc(e_incy); axis equal; axis tight;


%%% Phase attachment program
img_size1=size(e_incx);
img_size2=size(e_incy);
img_bin1=e_incx>0.7; %binarize
img_bin2=e_incy>0.7;
Img_bin1=e_incx>0.7; %binarize
Img_bin2=e_incy>0.7;

figure(1)
subplot(1,2,1);imagesc(e_incx);axis image;colormap jet; 
subplot(1,2,2);imagesc(e_incy);axis image;colormap jet; 
axis xy
hold on

%% Linear fitting
figure(2)
col1 = [1986.*ratio 2013.*ratio];
row1 = [2027.*ratio 1970.*ratio];
col2 = [1975.*ratio 2027.*ratio];
row2 = [1993.*ratio 2008.*ratio];



% Calculation of fitting function
degree = 1; % fitting order
func_parameter1 = polyfit(col1, row1, degree);
func_parameter1(1)=func_parameter1(1).*til1;
func_parameter1(2)=func_parameter1(2)+itc1;
slope1=func_parameter1(1);
intercept1=func_parameter1(2);

func_parameter2 = polyfit(col2, row2, degree);
func_parameter2(1)=func_parameter2(1).*til2;
func_parameter2(2)=func_parameter2(2)+itc2;
slope2=func_parameter2(1);
intercept2=func_parameter2(2);

% 近似直線のx座標を生成
xFit1 = 0:img_size1(1);
xFit2 = 0:img_size2(1);
% 近似直線のy座標を計算
yFit1 = polyval(func_parameter1, xFit1);
yFit2 = polyval(func_parameter2, xFit2);

[row1, col1] = find(img_bin1);
[row2, col2] = find(img_bin2);
subplot(1,2,1);scatter(col1,row1,'.');
hold on
plot(xFit1, yFit1, 'b-', 'LineWidth', 2);
axis([0 L 0 L]);
axis square;

subplot(1,2,2);scatter(col2,row2,'.');
hold on
plot(xFit2, yFit2, 'b-', 'LineWidth', 2);
axis([0 L 0 L]);
axis square;




%% Generation of phase boundary
figure(3)
% Intensity location
for i=1:n1+1
    I=i-1;
    interval1=(pointx_n-pointx_0)/n1;
    point_x1(i)=pointx_0+I*interval1;
    point_y1(i)=point_x1(i)*slope1+intercept1;
end
for i=1:n2+1
    I=i-1;
    interval2=(pointy_n-pointy_0)/n2;
    point_x2(i)=pointy_0+I*interval2;
    point_y2(i)=point_x2(i)*slope2+intercept2;
end

% Midpoint of intensity interval
boundary_x1(1)=2000.*ratio; boundary_y1(1)=boundary_x1(1)*slope1+intercept1;
boundary_x2(1)=2000.*ratio; boundary_y2(1)=boundary_x2(1)*slope2+intercept2;

%boundary_x1(2)=2000.*ratio; boundary_y1(2)=boundary_x1(2)*slope1+intercept1;
%boundary_x2(2)=2673.*ratio; boundary_y2(2)=boundary_x2(2)*slope2+intercept2;

%boundary_x1(3)=2714; boundary_y1(3)=boundary_x1(3)*slope1+intercept1;
%boundary_x2(3)=2722.*ratio; boundary_y2(3)=boundary_x2(3)*slope2+intercept2;

% boundary_x1(4)=350; boundary_y1(4)=boundary_x1(4)*slope1+intercept1;
%boundary_x2(4)=2770.*ratio; boundary_y2(4)=boundary_x2(4)*slope2+intercept2;

% Phase boundary
slope_inv1=-1/slope1;
slope_inv2=-1/slope2;

subplot(1,2,1);imagesc(Img_bin1);axis image;colormap jet; 
hold on
scatter(boundary_x1,boundary_y1);
for i=1:n1
    bound1=slope_inv1.*xFit1-slope_inv1.*boundary_x1(i)+boundary_y1(i);
    plot(xFit1,bound1)
end
axis xy;axis image

subplot(1,2,2);imagesc(Img_bin2);axis image;colormap jet; 
hold on
scatter(boundary_x2,boundary_y2);
for i=1:n2
    bound2=slope_inv2.*xFit1-slope_inv2.*boundary_x2(i)+boundary_y2(i);
    plot(xFit2,bound2)
end
axis xy;axis image

%% phase plot
phase1=ones(img_size1(1)).*exp(1j.*2.*pi);
phase2=ones(img_size2(1)).*exp(1j.*2.*pi);
for i=1:n1+1
    if i==1
        for x=1:img_size1(1)
            for y=1:img_size1(1)
                if y>slope_inv1*x-slope_inv1*boundary_x1(1)+boundary_y1(1)
                    phase1(y,x)=phase1(y,x).*exp(1j*pi);
                end
            end
        end

    elseif i<n1+1
        if 1<i
            if rem(i,2)==1
                for x=1:img_size1(1)
                    for y=1:img_size1(1)
                        if slope_inv1*x-slope_inv1*boundary_x1(i)+boundary_y1(i)<y
                            if y<slope_inv1*x-slope_inv1*boundary_x1(i-1)+boundary_y1(i-1)
                                phase1(y,x)=phase1(y,x).*exp(1j*pi);
                            end
                        end
                    end
                end
            end
        end

    else
        if rem(i,2)==1
            for x=1:img_size1(1)
                for y=1:img_size1(1)
                    if slope_inv1*x-slope_inv1*boundary_x1(n1)+boundary_y1(n1)>y
                        phase1(y,x)=phase1(y,x).*exp(1j*pi);
                    end
                end
            end
        end
     end
end


for i=1:n2+1
    if i==1
        for x=1:img_size2(1)
            for y=1:img_size2(1)
                if y<slope_inv2*x-slope_inv2*boundary_x2(1)+boundary_y2(1)
                    phase2(y,x)=phase2(y,x).*exp(1j*pi);
                end
            end
        end

    elseif i<n2+1
        if 1<i
            if rem(i,2)==1
                for x=1:img_size2(1)
                    for y=1:img_size2(1)
                        if slope_inv2*x-slope_inv2*boundary_x2(i)+boundary_y2(i)>y
                            if y>slope_inv2*x-slope_inv2*boundary_x2(i-1)+boundary_y2(i-1)
                                phase2(y,x)=phase2(y,x).*exp(1j*pi);
                            end
                        end
                    end
                end
            end
        end

    else
        if rem(i,2)==1
            for x=1:img_size2(1)
                for y=1:img_size2(1)
                    if slope_inv2*x-slope_inv2*boundary_x2(n2)+boundary_y2(n2)<y
                        phase2(y,x)=phase2(y,x).*exp(1j*pi);
                    end
                end
            end
        end
     end
end

        
figure(4)
subplot(1,2,1);imagesc(angle(phase1));axis image;axis xy;colormap jet; 
subplot(1,2,2);imagesc(angle(phase2));axis image;axis xy;colormap jet; 
axis xy







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%initial parameters
%L=2^11; % main domain size2048

z0=0; %interface position (focus at z=0)
z=linspace(-2,2,41); %axial domain
m0=1/2; %coordinate shift
xmax=10; %half the size of the cropped output in wavelengths
%(final size of images is 2xmax)
n1=1; % refractive index for incident field (immersion oil)
n2=1; % refractive index for reflected field (immersion oil)
n3=1; % refractive index for transmitted field (after interface)
e1=n2^2;
e2=n3^2;
m1=1; %permeability
m2=1; %permeability
f=1800; % effective focal length of the lens in um (10ˆ−6 m)
NA=0.99; %effective numerical aperture
lam=0.64; %wavelength in um
k=2*pi/lam; %Wave number
l=0.64; %wavelength um
R=f*NA/n2; %aperture radius in um;
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
k1=n2*k0;
k2=n3*k0;
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

PhaseShift=m0*2*pi*X/L+m0*2*pi*Y/L+0.8; %correction phase
%PhaseShift=; %correction phase
PhaseShift=PhaseShift.*0;
% figure(1)
% imagesc(PhaseShift)
%Center of Fourier transform displaced (vert and hor) by m0.
%shift of \∆ kx=\∆ ky= m0*(2*pi/L).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input intensity
%%%%%%%%%%%%%%%%%%%%%input field initial field transverse


% p=3;m=3;ec=2;lam=1;N=L+1;
% [IGB,Xin,Yin]=Ince_Gaussian(15e-3,N,0,p,m,ec,2e-4,(2*pi/lam),0);
% e_incy=rotate(IGB,90);
% e_incy=e_incy(1:L,1:L);

phi_incx=angle(phase1);
phi_incy=angle(phase2);
% E_incx=abs(E_incx.*conj(E_incx)); E_incy=abs(E_incy.*conj(E_incy));

E_incx=e_incx.^2./max(max(e_incx.^2)); E_incy=e_incy.^2./max(max(e_incy.^2));



figure(5)
fontsize(5,25,"points");colormap('jet');
subplot(2,2,1);imagesc(E_incx);axis image; colorbar(); title('Input intensity(x)')
subplot(2,2,2);imagesc(E_incy);axis image; colorbar(); title('Input intensity(y)')
subplot(2,2,3);imagesc(phi_incx);axis image;colorbar(); clim([-pi pi]); title('Input phase(x)') 
subplot(2,2,4);imagesc(phi_incy);axis image;colorbar(); clim([-pi pi]); title('Input phase(y)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%round soft aperture (spatial coordinates)
Theta=0.5*(1+tanh((1.5/1)*(N-sqrt((x_inf/dx_inf).^2+(y_inf/dx_inf).^2))));


%angular space
% kmax=1.5.*R*k1/f;
% Theta=0.5*(1+tanh((1.5/dk)*(kmax-sqrt((kx).^2+(ky).^2))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%polarization
%cx=ones(L,L); cy=zeros(L,L);   % lin horizontal
%cx=zeros(L,L); cy=ones(L,L); % lin vertical
%cx=1*ones(L,L); cy=1*ones(L,L); % lin diag
cx=ones(L,L); cy=1i*ones(L,L); %%right circ
%cx=1i*ones(L,L); cy=ones(L,L); %%left circ
%cx=(kx)./sqrt(kx.ˆ2+ky.ˆ2); cy=(ky)./sqrt(kx.ˆ2+ky.ˆ2); %%radial cx=cos(phi) cy=sin(phi)
%cx=−ky./sqrt(kx.ˆ2+ky.ˆ2); cy= kx./sqrt(kx.ˆ2+ky.ˆ2);%%azimuthal cx=cos(phi) cy=sin(phi)
%phik=atan2(ky,kx); s12=8; cx=cos((s12/2)*phik); cy=sin((s12/2)*phik);%flower
%phik=atan2(ky,kx); s12=−8; cx=cos((s12/2)*phik); cy=sin((s12/2)*phik);% spider

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Factors for E inf
CF1=-(1./kz1).*sqrt(kz1./kM1).*(1./(kx.^2+ky.^2))*1i*f*exp(-1i*f*2*pi*n2/lam); %Factors for E_inf_r and f
CF2=-(1./kz2).*sqrt(kz1./kM1).*(1./(kx.^2+ky.^2))*1i*f*exp(-1i*f*2*pi*n2/lam).*(kz2./kz1); %Factors for E inf t

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
% figure(2)
% imagesc(abs(Hzf))

%%%%%%%%%% focal
crop=L/25;
%x component E inf
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
Ex_n=Ex./max(max(Et));
Ey_n=Ey./max(max(Et));
Ez_n=Ez./max(max(Et));
Et_n=Et./max(max(Et));
%%%%%%%%%%%

figure(6)
subplot(1,3,1);imagesc(abs(Fieldxxf.*conj(Fieldxxf)));axis image; colorbar(); title('Total Intensity');
subplot(1,3,2);imagesc(abs(Fieldxyf.*conj(Fieldxyf)));axis image; colorbar(); title('x')
subplot(1,3,3);imagesc(abs(Fieldxzf));axis image; colorbar(); title('y')
round(max(max(Ex_n)),2)
figure(7)
subplot(1,4,1);imagesc(Et_n);axis image;  title('Total Intensity');axis off; clim([0 1]); %colorbar('Ticks',[0 1],'TickLabels',{0,1});
subplot(1,4,2);imagesc(Ex_n);axis image;  title('x');axis off; clim([0 round(max(max(Ex_n)),2)]); %colorbar('Ticks',[0 round(max(max(Ex_n)),2)],'TickLabels',{0,round(max(max(Ex_n)),2)});
subplot(1,4,3);imagesc(Ey_n);axis image;  title('y');axis off; clim([0 round(max(max(Ey_n)),3)]); %colorbar('Ticks',[0 round(max(max(Ey_n)),3)],'TickLabels',{0,round(max(max(Ey_n)),3)});
subplot(1,4,4);imagesc(Ez_n);axis image;  title('z');axis off; clim([0 round(max(max(Ez_n)),2)]); %colorbar('Ticks',[0 round(max(max(Ez_n)),2)],'TickLabels',{0,round(max(max(Ez_n)),2)});
fontsize(7,22,"points");colormap('hot');
figure(8)
subplot(1,4,2);imagesc(angle(Exf));axis image; colorbar(); clim([-pi pi]); title('phase x');axis off;
subplot(1,4,3);imagesc(angle(Eyf));axis image; colorbar(); clim([-pi pi]); title('phase y');axis off;
subplot(1,4,4);imagesc(angle(Ezf));axis image; colorbar(); clim([-pi pi]); title('phase z');axis off;
fontsize(8,22,"points");colormap('jet');
%%plot yz cross section through the center at xmax2+1
aa=angle(Exf);

%% 
function y=HGmode(n,m,x,y,w,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*exp(-((x.^2+y.^2)./W^2)-(1j*k.*(x.^2+y.^2)./(2*R))).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y);%.*exp(gouy);
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

function y=g_filter(img,Size,sigma)
    % ガウシアンフィルタのサイズと標準偏差を指定
    filterSize = Size; % フィルタサイズを調整してください
    sig = sigma;    % ガウシアン分布の標準偏差を調整してください
    
    % ガウシアンカーネルの作成
    [X, Y] = meshgrid(-floor(filterSize/2):floor(filterSize/2), -floor(filterSize/2):floor(filterSize/2));
    gaussianKernel = exp(-(X.^2 + Y.^2) / (2 * sig^2));
    gaussianKernel = gaussianKernel / sum(gaussianKernel(:));
    
    % 画像サイズの取得
    [rows, cols] = size(img);
    
    % 処理結果を格納するための配列
    filteredImage = zeros(size(img));
    
    % ガウシアンフィルタの適用
    for i = floor(filterSize/2)+1:rows-floor(filterSize/2)
        for j = floor(filterSize/2)+1:cols-floor(filterSize/2)
            % フィルターサイズの範囲内の画素値を取得
            neighborhood = img((i-floor(filterSize/2)):(i+floor(filterSize/2)), ...
                                    j-floor(filterSize/2):j+floor(filterSize/2));
            
            % ガウシアンカーネルと画素値の畳み込み
            filteredImage(i,j) = sum(neighborhood(:) .* gaussianKernel(:));
        end
    end
    y=filteredImage;
end

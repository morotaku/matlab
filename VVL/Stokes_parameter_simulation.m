close all;clc;clear;

%%%% s1 import data to matlab from folder
%%% x-y　coordinate
N=501;
L=3;
%L=35; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi,r] = cart2pol(x,y);

%%% Parameter
w=1;
lam=1;
NA=0.9;
f=5*10^3;
w0=0.61*lam/NA;
zr=pi*w^2/lam;
k=2*pi/lam;
z=zr.*0;
alpha=asin(NA);
thitaP=phi;
R=z+zr^2./z;
n1=4;n2=3;
p=3;m=3;ec=3;
[IGB,Xin,Yin]=Ince_Gaussian(15e-3,501,0,p,m,ec,5e-3,(2*pi/lam),0);
% IGB=abs(IGB.^2);
e_h=hgmode(0,4,x,y,w,z,R,k);

i_h=abs(e_h.*conj(e_h));
% e_v=IGB;
% e_v=rotate(e_v,45);
e_v=hgmode(0,3,x,y,w,z,R,k)+8.*hgmode(2,1,x,y,w,z,R,k)+8.*hgmode(1,2,x,y,w,z,R,k)+hgmode(3,0,x,y,w,z,R,k);
i_v=abs(e_v.*conj(e_v));
I00 = i_h./max(max(i_h));
I90 = i_v./max(max(i_v));
I45 = I00+I90;
I135 = I00+I90;


I00 = I00./max(max(I00));
I90 = I90./max(max(I90));
I45 = I45./max(max(I45));
I135 = I135./max(max(I135));



%%%% s2 range of image 0 to 255
range=255;
mx=max(max(I00+I90));

% I00=(I00./mx)*range;
% I90=(I90./mx)*range;
% I45=(I45./mx)*range;
% I135=(I135./mx)*range;


%%%%%%%  s3 stokes parameter
s0=(I00+I90);
s1=(I00-I90);%./(I00+I90);
s2=(I45-I135);%./(I135+I45);
s3=sqrt(s0.^2-s1.^2-s2.^2);

s0_n=s0./max(max(s0));
s1_n=s1./max(max(s0));
s2_n=s2./(max(max(s0)));
s3_n=s3./max(max(s0));

phase1=angle(e_h);
phase2=angle(e_v);
subplot(1,2,1);imagesc(phase1);colorbar();axis xy;
subplot(1,2,2);imagesc(phase2);colorbar();axis xy;
phase2=angle(-e_v.*exp(1j*pi/2));

%imagesc(phase1);colorbar();axis xy
% subplot(1,2,1);imagesc(phase1);axis image;axis xy;colormap jet; colorbar;clim([-1.5 3])
% subplot(1,2,2);imagesc(phase2);axis image;axis xy;colormap jet; colorbar;clim([-1.5 3])
phase_dif=phase2-phase1;

right_cir=sin(phase_dif)>0;
left_cir=sin(phase_dif)<0;
s3_sign=right_cir-left_cir;

s3_n=s3_n.*s3_sign;
%%%%%%  s4 pol angles
shi=0.5*atan2(s2,s1); % -90 to 90
xi=0.5*atan2(s3,(sqrt((s1.^2)+ (s2.^2)))); % -45 to 45
ximod=[-pi/2*ones(size(xi,1),1),xi,pi/2*ones(size(xi,1),1)];


%figure
figure(1)
imagesc(s0_n);axis image; axis off;colormap jet; title('s0');
colorbar();clim([-1 1]);

figure(2)
imagesc(s1_n);axis image; axis off;colormap jet; title('s1');
colorbar();clim([-1 1]);

figure(3)
imagesc(s2_n);axis image; axis off;colormap jet; title('s2');
colorbar();clim([-1 1]);

figure(4)
imagesc(s3_n); axis image;axis off;colormap jet;  title('s3');
colorbar();clim([-1 1]);

% %figure
% figure(1)
% subplot(1,4,1);imagesc(s0);axis image; axis off;colormap jet; title('s0');
% colorbar()
% figure(2)
% subplot(1,4,2);imagesc(s1);axis image; axis off; title('s1');
% colorbar()
% figure(3)
% subplot(1,4,3);imagesc(s2);axis image; axis off; title('s2');
% colorbar()
% figure(4)
% subplot(1,4,4);imagesc(s3); axis image;axis off;  title('s3');
% colorbar()
% saveas(gcf,'S1-3_re.png')


figure(2)
subplot(1,2,1);imagesc(shi);axis image; axis off;colormap jet; title('shi');
subplot(1,2,2);imagesc(ximod);axis image; axis off; title('xi');
saveas(gcf,'shi and xi_re.png')
figure(3)
f_stokes_plot(s1,s2,s3,50);

function y=HGmode(n,m,x,y,w,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-((x.^2+y.^2)./W^2)-(1j*k.*(x.^2+y.^2)./R));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=U_n;%./max(max(abs(u)));
end
function y=hgmode(n,m,x,y,w,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-((x.^2+y.^2)./W^2)-(1j*k.*(x.^2+y.^2)./R));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));
end

function y=HermitePol(n,x)
Hn=0;
for m=0:floor(n/2)
    Hn=Hn+((-1)^m.*(2.*x).^(n-2*m))./(factorial(m)*factorial(n-2*m));
end
Hn=Hn.*factorial(n);
y=Hn;
end

function u=HGmode45(n,x,y,r,w,z,zr,k,R)
    HG=0;
    for i=0:n
        HG=HG+nchoosek(n,i).*hgmode(i,n-i,x,y,w,z,R,k);
    end
    u=HG.*sqrt(2)^n;
end

function u=HGmode_45(n,x,y,r,w,z,zr,k,R)
    HG=0;
    for i=0:n
        if rem(i,2)==0
            HG=HG-nchoosek(n,i).*hgmode(i,n-i,x,y,w,z,R,k);
        else
            HG=HG+nchoosek(n,i).*hgmode(i,n-i,x,y,w,z,R,k);
        end
    end
    u=HG.*sqrt(2)^n;
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


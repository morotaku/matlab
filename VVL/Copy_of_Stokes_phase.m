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


e_h=hgmode(0,3,x,y,w,z,R,k);
i_h=abs(e_h.*conj(e_h));
e_v=hgmode(3,0,x,y,w,z,R,k);
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
s0=(abs(e_h).^2+abs(e_v).^2);
s1=(abs(e_h).^2-abs(e_v).^2);
s2=(conj(e_v).*e_h+conj(e_h).*e_v);
s3=1j.*(conj(e_v).*e_h-conj(e_h).*e_v);

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

% Stp12=atan(s2_n./s1_n);
% Stp23=atan(s3_n./s2_n);
% Stp31=atan(s1_n./s3_n);
Stp12=atan2(s2_n,s1_n);
Stp23=atan2(s3_n,s2_n);
Stp31=atan2(s1_n,s3_n);

Stp12(Stp12<0)=Stp12(Stp12<0)+2.*pi;
Stp23(Stp23<0)=Stp23(Stp23<0)+2.*pi;
Stp31(Stp31<0)=Stp31(Stp31<0)+2.*pi;



%figure
figure(1)
subplot(1,4,1);imagesc(s0_n);axis image; axis off;colormap jet; title('s0');colorbar();clim([-1 1]);
subplot(1,4,2);imagesc(s1_n);axis image; axis off;colormap jet; title('s1');colorbar();clim([-1 1]);
subplot(1,4,3);imagesc(s2_n);axis image; axis off;colormap jet; title('s2');colorbar();clim([-1 1]);
subplot(1,4,4);imagesc(s3_n); axis image;axis off;colormap jet;  title('s3');colorbar();clim([-1 1]);

figure(2)
colormap('hsv');
subplot(1,3,1);imagesc(Stp12);axis image; axis off;colormap jet; title('Φ12');colorbar();clim([0 2*pi]);
subplot(1,3,2);imagesc(Stp23);axis image; axis off;colormap jet; title('Φ23');colorbar();clim([0 2*pi]);
subplot(1,3,3);imagesc(Stp31);axis image; axis off;colormap jet; title('Φ31');colorbar();clim([0 2*pi]);
% subplot(1,3,1);imagesc(Stp12);axis image; axis off;colormap jet; title('Φ12');colorbar();clim([-pi pi]);
% subplot(1,3,2);imagesc(Stp23);axis image; axis off;colormap jet; title('Φ23');colorbar();clim([-pi pi]);
% subplot(1,3,3);imagesc(Stp31);axis image; axis off;colormap jet; title('Φ31');colorbar();clim([-pi pi]);
fontsize(2,25,"points")

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
y=LG1./max(max(abs(LG1)));

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


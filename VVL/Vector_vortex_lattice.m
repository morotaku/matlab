clear all
close all
clc
L=2.25; %Display lange
p=3;m=3;ec=2;
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=L/8.5;
%%% Polarizer angle 
ang=input('polarizer angle(thita=-1 → no polarizer): '); 

%%% Topological charge
n1=0; %North pole
n2=1; %South pole
%%% Radial index
m1=2;
m2=0;


%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude

w=1; %Beam waist
lam=1; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size


%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;

X=linspace(-L,L,N);
Y=linspace(-L,L,N);

[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
hg1=hgmode(n1,m1,x,y,w,r,z,R,k);
HG1=abs(hg1.*conj(hg1));
HG1_n=HG1/max(max(HG1));
hg2=hgmode(n2,m2,x,y,w,r,z,R,k);
%[hg2,Xin,Yin]=Ince_Gaussian(15e-3,1001,0,p,m,ec,6e-3,(2*pi/lam),0);
%hg2=rotate(hg2,-45);
HG2=abs(hg2.*conj(hg2));
HG2_n=HG2/max(max(HG2));

figure(1)
colormap('hot')
subplot(1,5,1);imagesc(HG1_n+HG2_n);axis equal; axis tight; axis off;
subplot(1,5,2);imagesc(HG1_n);axis equal; axis tight; axis off;
subplot(1,5,3);imagesc(HG1_n+HG2_n);axis equal; axis tight; axis off;
subplot(1,5,4);imagesc(HG2_n);axis equal; axis tight; axis off;
subplot(1,5,5);imagesc(HG1_n+HG2_n);axis equal; axis tight; axis off;
axis equal; axis off; axis tight;
hold on


function y=HGmode_45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        if rem(i,2)==1
            hg=hg-nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        else
            hg=hg+nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        end
    end
    y=hg./(sqrt(2)^n);
    %y=hg./max(max(abs(hg)));
end

function y=HGmode45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        hg=hg+nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
    end
    y=hg./(sqrt(2)^n);
end

function y=hgmode_45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        if rem(i,2)==1
            hg=hg-nchoosek(n,i).*hgmode(i,n-i,x,y,w,r,z,R,k);
        else
            hg=hg+nchoosek(n,i).*hgmode(i,n-i,x,y,w,r,z,R,k);
        end
    end
    y=hg./(sqrt(2)^n);
    %y=hg./max(max(abs(hg)));
end

function y=hgmode45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        hg=hg+nchoosek(n,i).*hgmode(i,n-i,x,y,w,r,z,R,k);
    end
    y=hg./(sqrt(2)^n);
end

function y=HGmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-(r.^2./W^2)-(1j*k.*r.^2/(2*R)));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));
end

function y=hgmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-((x.^2+y.^2)./W^2));%-(1j*k.*(x.^2+y.^2)./(2*R))).*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));%./max(max(abs(u)));
end

function y=HermitePol(n,x)
Hn=0;
for m=0:floor(n/2)
    Hn=Hn+((-1)^m.*(2.*x).^(n-2*m))./(factorial(m)*factorial(n-2*m));
end
Hn=Hn.*factorial(n);
y=Hn;
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
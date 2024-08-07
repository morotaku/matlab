%%% Simple computer program to calculate arbitrary tightly focused (propagating andevanescent) vector light fields
clear all; close all; clc;
fig_size=80;zrange=1;
Linew=1;
scale=0.2;

Exf=importdata('tightly_focused_image/num3/Exf.mat'); Exf=imresize(Exf,[fig_size fig_size]);
Eyf=importdata('tightly_focused_image/num3/Eyf.mat'); Eyf=imresize(Eyf,[fig_size fig_size]);
Ezf=importdata('tightly_focused_image/num3/Ezf.mat'); Ezf=imresize(Ezf,[fig_size fig_size]);

Ex=abs(Exf.*conj(Exf)); % total field at the material n 2 incident+reflected
Ey=abs(Eyf.*conj(Eyf));
Ez=abs(Ezf.*conj(Ezf));
Et=Ex+Ey+Ez;
et=max(max(Et));
Ex_n=Ex./max(max(et));
Ey_n=Ey./max(max(et));
Ez_n=Ez./max(max(et));
Et_n=Et./max(max(et));

phase_x=angle(Exf);phase_y=angle(Eyf);phase_z=angle(Ezf);
%phase_x(phase_x>2.8)=-phase_x(phase_x>2.8); phase_y(phase_y>2.8)=-phase_y(phase_y>2.8); phase_z(phase_z>2.8)=-phase_z(phase_z>2.8);



%%%%%%%%%%%
% Polarization strip
Ran=size(Exf,1);
x2=linspace(0,Ran,Ran); %%%auxiliary coordinate for correction
y2=x2; [X2,Y2]=meshgrid(x2,y2);
[phi,r] = cart2pol(X2,Y2);

ex=real(Exf)+1j.*imag(Exf);
ey=real(Eyf)+1j.*imag(Eyf);
ez=real(Ezf)+1j.*imag(Ezf);
EE=dot(ex,ex)+dot(ey,ey)+dot(ez,ez);

% x-yPlane Projection
ax=real(ex.*sqrt(EE))./abs(sqrt(EE));
ax_n=ax./max(max(ax));
ay=real(ey.*sqrt(EE))./abs(sqrt(EE));
ay_n=ay./max(max(ay));
az=real(ez.*sqrt(EE))./abs(sqrt(EE));
az_n=az./max(max(az));


cnt1=0;
cnt2=0;
cnt3=0;
g1=[];
g2=[];
g3=[];
strip1=Ran/8.5;
strip2=Ran/4;
strip3=Ran/2.5;
for i=1:Ran
    for j=1:Ran
        I=i-(Ran+1)/2;
        J=j-(Ran+1)/2;
        if strip1*0.95<sqrt(I^2+J^2)
            if sqrt(I^2+J^2)<strip1*1.05
                cnt1=cnt1+1;
                g1(cnt1)=i+j*Ran;
                aX1(cnt1)=ax(i,j); aY1(cnt1)=ay(i,j); aZ1(cnt1)=az(i,j);
                x_1(cnt1)=X2(i,j); y_1(cnt1)=Y2(i,j); z_1(cnt1)=0;
            end
        end

        if strip2*0.975<sqrt(I^2+J^2)
            if sqrt(I^2+J^2)<strip2*1.025
                cnt2=cnt2+1;
                g2(cnt2)=i+j*Ran;
                aX2(cnt2)=ax(i,j); aY2(cnt2)=ay(i,j); aZ2(cnt2)=az(i,j);
                x_2(cnt2)=X2(i,j); y_2(cnt2)=Y2(i,j); z_2(cnt2)=0;
            end
        end

        if strip3*0.9825<sqrt(I^2+J^2)
            if sqrt(I^2+J^2)<strip3*1.0175
                cnt3=cnt3+1;
                g3(cnt3)=i+j*Ran;
                aX3(cnt3)=ax(i,j); aY3(cnt3)=ay(i,j); aZ3(cnt3)=az(i,j);
                x_3(cnt3)=X2(i,j); y_3(cnt3)=Y2(i,j); z_3(cnt3)=0;
            end
        end
    end
end


figure(1)
colormap('hot')
subplot(1,4,1);imagesc(Et_n);axis image; colorbar; clim([0 1]); title('Total Intensity'); axis off;
hold on; scatter(x_1,y_1,100,'white','.'); scatter(x_2,y_2,100,'blue','.'); scatter(x_3,y_3,100,'red','.');
subplot(1,4,2);imagesc(Ex_n);axis image; colorbar;axis off; clim([0 max(max(Ex_n))]); title('x');
hold on; scatter(x_1,y_1,100,'white','.'); scatter(x_2,y_2,100,'blue','.'); scatter(x_3,y_3,100,'red','.');
subplot(1,4,3);imagesc(Ey_n);axis image; colorbar;axis off; clim([0 max(max(Ey_n))]); title('y');
hold on; scatter(x_1,y_1,100,'white','.'); scatter(x_2,y_2,100,'blue','.'); scatter(x_3,y_3,100,'red','.');
subplot(1,4,4);imagesc(Ez_n);axis image; colorbar;axis off; clim([0 max(max(Ez_n))]); title('z');
hold on; scatter(x_1,y_1,100,'white','.'); scatter(x_2,y_2,100,'blue','.'); scatter(x_3,y_3,100,'red','.');
fontsize(1,22,"points");

figure(2)
colormap('hsv')
subplot(1,4,2);imagesc(phase_x);axis image; colorbar(); axis off; clim([-pi pi]);
hold on; scatter(x_1,y_1,100,'white','.'); scatter(x_2,y_2,100,'blue','.'); scatter(x_3,y_3,100,'red','.');
subplot(1,4,3);imagesc(phase_y);axis image; colorbar(); axis off; clim([-pi pi]);
hold on; scatter(x_1,y_1,100,'white','.'); scatter(x_2,y_2,100,'blue','.'); scatter(x_3,y_3,100,'red','.');
subplot(1,4,4);imagesc(phase_z);axis image; colorbar(); axis off; clim([-pi pi]);
hold on; scatter(x_1,y_1,100,'white','.'); scatter(x_2,y_2,100,'blue','.'); scatter(x_3,y_3,100,'red','.');
fontsize(2,22,"points");


%%%%%%% Polarization strip plot


%%% Plarization strip1
figure(3);

xrange1=strip1.*1.4; yrange1=strip1.*1.4; 
q1_p=quiver3(x_1,y_1,z_1,aX1.*scale,aY1.*scale,aZ1.*scale./2,1);
q1_p.LineWidth=Linew; q1_p.Color=[1 0 0];
hold on
q1_n=quiver3(x_1,y_1,z_1,-aX1.*scale,-aY1.*scale,-aZ1.*scale./2,1);
q1_n.LineWidth=Linew; q1_n.Color=[0 0 1];

% xyplane projection
z_plane=ones(1,numel(x_1)).*-zrange;
qxy1_p=quiver3(x_1,y_1,z_plane,aX1./5,aY1./5,aZ1.*0);
qxy1_p.LineWidth=Linew;qxy1_p.Color=[1 0.5 0.5];
hold on
qxy1_n=quiver3(x_1,y_1,z_plane,-aX1./5,-aY1./5,aZ1.*0,'b');
qxy1_n.LineWidth=Linew;qxy1_n.Color=[0.5 0.5 1];

xlim([Ran/2-xrange1 Ran/2+xrange1]);ylim([Ran/2-yrange1 Ran/2+yrange1]);zlim([-zrange zrange]);
xlabel('x');ylabel('y');zlabel('z');
title('Polarization strip 1');fontsize(3,18,"points");
axis off;

%%% Plarization strip2
figure(4)

xrange2=strip2*1.4; yrange2=strip2.*1.4;
q2_p=quiver3(x_2,y_2,z_2,aX2.*scale,aY2.*scale,aZ2.*scale./2,1);
q2_p.LineWidth=Linew; q2_p.Color=[1 0 0];
hold on
q2_n=quiver3(x_2,y_2,z_2,-aX2.*scale,-aY2.*scale,-aZ2.*scale./2,1);
q2_n.LineWidth=Linew; q2_n.Color=[0 0 1];

% xyplane projection
z_plane=ones(1,numel(x_2)).*-zrange;
qxy2_p=quiver3(x_2,y_2,z_plane,aX2./5,aY2./5,aZ2.*0);
qxy2_p.LineWidth=Linew;qxy2_p.Color=[1 0.5 0.5];
hold on
qxy2_n=quiver3(x_2,y_2,z_plane,-aX2./5,-aY2./5,aZ2.*0,'b');
qxy2_n.LineWidth=Linew;qxy2_n.Color=[0.5 0.5 1];

xlim([Ran/2-xrange2 Ran/2+xrange2]);ylim([Ran/2-yrange2 Ran/2+yrange2]);zlim([-zrange zrange]);
xlabel('x');ylabel('y');zlabel('z');
title('Polarization strip 2');fontsize(4,18,"points");
axis off;

%%% Plarization strip3
figure(5)

xrange3=strip3*1.4; yrange3=strip3*1.4; 
q3_p=quiver3(x_3,y_3,z_3,aX3.*scale,aY3.*scale,aZ3.*scale./2,1);
q3_p.LineWidth=Linew; q3_p.Color=[1 0 0];
hold on
q3_n=quiver3(x_3,y_3,z_3,-aX3.*scale,-aY3.*scale,-aZ3.*scale./2,1);
q3_n.LineWidth=Linew; q3_n.Color=[0 0 1];

% xyplane projection
z_plane=ones(1,numel(x_3)).*-zrange;
qxy3_p=quiver3(x_3,y_3,z_plane,aX3./5,aY3./5,aZ3.*0);
qxy3_p.LineWidth=Linew;qxy3_p.Color=[1 0.5 0.5];
hold on
qxy3_n=quiver3(x_3,y_3,z_plane,-aX3./5,-aY3./5,aZ3.*0,'b');
qxy3_n.LineWidth=Linew;qxy3_n.Color=[0.5 0.5 1];

xlim([Ran/2-xrange3 Ran/2+xrange3]);ylim([Ran/2-yrange3 Ran/2+yrange3]);zlim([-zrange zrange]);
xlabel('x');ylabel('y');zlabel('z');
title('Polarization strip 3');fontsize(5,18,"points");
axis off;


%% Function
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
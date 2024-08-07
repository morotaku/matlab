clear all
close all
clc
n1=4;
n2=3;
%%%傾き補正
til1=1.75;
til2=1;
%%%切片補正
itc1=230;
itc2=0;
pointx_0=205; %左端のローブ
pointx_n=405;
pointy_0=135; %左端のローブ
pointy_n=505;



%% import data and binarize
img1=importdata('6/0.csv');
img2=importdata('6/90.csv');
img1=img1./max(max(img1));
img2=img2./max(max(img2));
img_size1=size(img1);
img_size2=size(img2);
img_bin1=img1>0.8; %binarize
img_bin2=img2>0.6;

figure(1)
subplot(1,2,1);imagesc(img1);axis image;colormap jet; grid on;colorbar;axis xy
subplot(1,2,2);imagesc(img2);axis image;colormap jet; grid on;colorbar;axis xy
hold on

%% Linear fitting
figure(2)
[row1, col1] = find(img_bin1);
[row2, col2] = find(img_bin2);


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

subplot(1,2,1);scatter(col1,row1,'.');
hold on
plot(xFit1, yFit1, 'b-', 'LineWidth', 2);
axis image;xlim([0 img_size1(1)]);ylim([0 img_size1(1)]);

subplot(1,2,2);scatter(col2,row2,'.');
hold on
plot(xFit2, yFit2, 'b-', 'LineWidth', 2);
axis image;xlim([0 img_size1(1)]);ylim([0 img_size1(1)]);




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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boundary_x1(1)=190; boundary_y1(1)=boundary_x1(1)*slope1+intercept1;
boundary_x2(1)=160; boundary_y2(1)=boundary_x2(1)*slope2+intercept2;

boundary_x1(2)=212.5; boundary_y1(2)=boundary_x1(2)*slope1+intercept1;
boundary_x2(2)=225; boundary_y2(2)=boundary_x2(2)*slope2+intercept2;
% 
boundary_x1(3)=232.5; boundary_y1(3)=boundary_x1(3)*slope1+intercept1;
boundary_x2(3)=300; boundary_y2(3)=boundary_x2(3)*slope2+intercept2;
% 
boundary_x1(4)=262.5; boundary_y1(4)=boundary_x1(4)*slope1+intercept1+15;
% boundary_x2(4)=337.5; boundary_y2(4)=boundary_x2(4)*slope2+intercept2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:n1
%     boundary_x1(i)=(point_x1(i)+point_x1(i+1))/2;
%     boundary_y1(i)=(point_y1(i)+point_y1(i+1))/2;
% end
% for i=1:n2
%     boundary_x2(i)=(point_x2(i)+point_x2(i+1))/2;
%     boundary_y2(i)=(point_y2(i)+point_y2(i+1))/2;
% end

% Phase boundary
slope_inv1=-1/slope1;
slope_inv2=-1/slope2;

subplot(1,2,1);imagesc(img1);axis image;colormap gray; axis xy;
xlim([0 img_size1(1)]);ylim([0 img_size1(1)]);
hold on
scatter(boundary_x1,boundary_y1);
for i=1:n1
    bound1=slope_inv1.*xFit1-slope_inv1.*boundary_x1(i)+boundary_y1(i);
    plot(xFit1,bound1)
end
axis xy;axis image;xlim([0 img_size1(1)]);ylim([0 img_size1(1)]);

subplot(1,2,2);imagesc(img2);axis image;colormap gray; axis xy;
xlim([0 img_size1(1)]);ylim([0 img_size1(1)]);
hold on
scatter(boundary_x2,boundary_y2);
for i=1:n2
    bound2=slope_inv2.*xFit1-slope_inv2.*boundary_x2(i)+boundary_y2(i);
    plot(xFit2,bound2)
end
axis xy;axis image;
xlim([0 img_size1(1)]);ylim([0 img_size1(1)]);

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

p=3;m=3;ec=2;lam=1;N=img_size2(1)+1;
[IGB,Xin,Yin]=Ince_Gaussian(15e-3,N,0,p,m,ec,4.5e-3,(2*pi/lam),0);
imagesc(abs(IGB));axis image;
IGB=rot90(IGB,-1);
aa=angle(IGB.*exp(1j*pi/2));%.*exp(-1j*17.*pi/180);
% figure(4)
% imagesc(aa);axis xy;
d=25;
for i=1:N
    for j=1:N-d
        aa(N+1-j,i)=aa(N+1-j-d,i);
    end
    aa(1:d,i)=aa(d+1,i);
end
aa(1,:)=[];
aa(:,1)=[];
% figure(5)
% imagesc(aa);axis xy;
% 回転する角度を指定 (例: 45度)
rotationAngle = 17;

% 画像のサイズを取得
[rows, cols, ~] = size(aa);

% 回転行列の作成
theta = deg2rad(rotationAngle);
rotationMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% 中心座標の計算
center = [(cols+1)/2, (rows+1)/2];

% 回転行列を用いて画像を回転
rotatedImage = zeros(size(aa));
for i = 1:rows
    for j = 1:cols
        % 回転前の座標を計算
        originalCoord = rotationMatrix * ([j; i] - center') + center';
        
        % 回転前の座標が画像の範囲内か確認してから値を代入
        if all(originalCoord > 0) && all(originalCoord <= [cols; rows])
            % バイリニア補間を行う場合は、適切な補間関数を実装してください
            
            rotatedImage(i, j, :) = aa(ceil(originalCoord(2)), ceil(originalCoord(1)), :);
        end
    end
end


figure(4)
subplot(1,2,1);imagesc(angle(phase1));axis image;axis xy;colormap jet; colorbar;
subplot(1,2,2);imagesc(rotatedImage);axis image;axis xy;colormap jet; colorbar;
axis xy

%% stokes parameter s3
figure(5)
%phase1=phase1.*exp(-1j*pi/2);
%imagesc(angle(phase1));colorbar; axis xy
phase_dif=(rotatedImage)-angle(phase1);
%phase_dif=g_filter(phase_dif,51,10);

pol45q=sin(phase_dif)+1;
pol135q=-sin(phase_dif)+1;



pol45q=pol45q./max(max(pol45q));
pol135q=pol135q./max(max(pol135q));

% subplot(1,2,1);imagesc(pol45q);colorbar;
% subplot(1,2,2);imagesc(pol135q);colorbar;
I00=img1;
I90=img2;
I45=importdata('6/45.csv');
I135=importdata('6/135.csv');

I00 = I00./max(max(I00));
I90 = I90./max(max(I90));
I45 = I45./max(max(I45));
I135 = I135./max(max(I135));

s0=(I00+I90);
s1=(I00-I90);%./(I00+I90);
s2=(I45-I135);%./(I135+I45);
s3=sqrt(abs(s0.^2-s1.^2-s2.^2));

s0_n=s0./max(max(s0));
s1_n=s1./max(max(s0));
s2_n=s2./(max(max(s0)));
s3_n=s3./max(max(s0));

right_cir=sin(phase_dif)>0;
left_cir=sin(phase_dif)<0;
s3_sign=right_cir-left_cir;

s3_sign=g_filter(s3_sign,51,10);
s3_n=s3_n.*s3_sign;
%s3_n=g_filter(s3_n,31,2);
imagesc(s3_n); axis image;axis off; colormap jet;
colorbar();clim([-1 1]);

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
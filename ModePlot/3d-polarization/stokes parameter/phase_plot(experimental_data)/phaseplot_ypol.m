clear all
close all
clc
n=4;
point_0=125; %左端のローブ
point_n=515;

%%%% s1 import data to matlab from folder
img=importdata('2.csv');
img=img./max(max(img));
img_size=size(img);
img_bin=img>0.475; %binarize

figure(1)
imagesc(img_bin)
axis xy
hold on
[row, col] = find(img_bin);
figure(2)
scatter(col,row,'.')
hold on
% 近似直線を計算
degree = 1; % フィッティング関数の次数
func_parameter = polyfit(col, row, degree);
slope=func_parameter(1);
intercept=func_parameter(2);
% 近似直線のx座標を生成
xFit = 0:img_size(1);
% 近似直線のy座標を計算
yFit = polyval(slope, xFit);
%intercept=polyval(slope, 0);
plot(xFit, yFit, 'b-', 'LineWidth', 2);
xlim([0 600])
ylim([0 600])



figure(3)

for i=1:n+1
    I=i-1;
    interval=(point_n-point_0)/n;
    point_x(i)=point_0+I*interval;
    point_y(i)=point_x(i)*slope+intercept;
end

for i=1:n
    boundary_x(i)=(point_x(i)+point_x(i+1))/2;
    boundary_y(i)=(point_y(i)+point_y(i+1))/2;
end

slope_inv=-1/slope;

imagesc(img_bin)
hold on
scatter(boundary_x,boundary_y);
for i=1:n
    bound=slope_inv.*xFit-slope_inv.*boundary_x(i)+boundary_y(i);
    plot(xFit,bound)
end
axis xy


phase=zeros(img_size(1));
for i=1:n+1
    if i==1
        for x=1:img_size(1)
            for y=1:img_size(1)
                if y<slope_inv*x-slope_inv*boundary_x(1)+boundary_y(1)
                    phase(y,x)=1j*pi;
                end
            end
        end

    elseif i<n+1
        if 1<i
            if rem(i,2)==1
                for x=1:img_size(1)
                    for y=1:img_size(1)
                        if slope_inv*x-slope_inv*boundary_x(i)+boundary_y(i)>y
                            if y>slope_inv*x-slope_inv*boundary_x(i-1)+boundary_y(i-1)
                                phase(y,x)=1j*pi;
                            end
                        end
                    end
                end
            end
        end

    else
        if rem(i,2)==1
            for x=1:img_size(1)
                for y=1:img_size(1)
                    if slope_inv*x-slope_inv*boundary_x(n)+boundary_y(n)<y
                        phase(y,x)=1j*pi;
                    end
                end
            end
        end
     end
end

        
figure(4)

imagesc(angle(phase))
colorbar
axis xy
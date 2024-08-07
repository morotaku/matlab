%%% Simple computer program to calculate arbitrary tightly focused (propagating andevanescent) vector light fields
clear all; close all; clc;
M=5;n1=4; n2=4;
N1=2800; % number of points for the aperture radius
N2=2800; % number of points for the aperture radius
range=110000; %display range
L=3000;
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


imp1=importdata('1/0.csv'); imp1_r=imresize(imp1,[fig_size fig_size]);
imp_x=abs(sqrt(imp1_r)); imp_x=imp_x./max(max(imp_x));
%imp_x=g_filter(imp_x,10,5);

imp2=importdata('1/90.csv'); imp2_r=imresize(imp2,[fig_size fig_size]);
imp2_r
imp_y=abs(sqrt(imp2_r)); imp_y=imp_y./max(max(imp_y));
%imp_y=g_filter(imp_y,10,5);
L2=size(imp_x);
L2=L2(1);
L3=5;
L4=-5;

adj=length(L/2-L2/2:L/2+L2/2)-L2;
e_incx=zeros(L); e_incy=zeros(L);
e_incx(L/2-L2/2+L4:L/2+L2/2-adj+L4,L/2-L2/2:L/2+L2/2-adj)=sqrt(imp_x);
% e_incy(L/2-L2/2:L/2+L2/2-adj,L/2-L2/2:L/2+L2/2-adj)=sqrt(imp_y);
e_incy(L/2-L2/2:L/2+L2/2-adj,L/2-L2/2-L3:L/2+L2/2-adj-L3)=sqrt(imp_y);
figure(1);imagesc(e_incx); axis equal; axis tight;
figure(2);imagesc(e_incy); axis equal; axis tight;


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
close all;clc;clear;

%%%% s1 import data to matlab from folder
I00 = importdata('lattice_mode/26 2_0001.ascii.csv');
I90 = importdata('lattice_mode/26 4_0001.ascii.csv');
I45 = importdata('lattice_mode/26 3_0001.ascii.csv');
I135 = importdata('lattice_mode/26 5_0001.ascii.csv');
test=importdata('lattice_mode/1.csv');
imagesc(test)
% I45q = importdata;%('I45q_0001.ascii.csv');
% I135q = importdata('I135q_0001.ascii.csv');
%%%% s2 range of image 0 to 255
range=255;
mx=max(max(I00+I90));

I00=(I00/mx)*range;
I90=(I90/mx)*range;
I45=(I45/mx)*range;
I135=(I135/mx)*range;
% I45q=(I45q/mx)*range;
% I135q=(I135q/mx)*range;

%%%%%%%  s3 stokes parameter
s0=(I00+I90);
s1=(I00-I90)./(I00+I90);
s2=(I45-I135)./(I135+I45);
% s3=(I45q-I135q)./(I135q+I45q);

%%%%%%  s4 pol angles
shi=0.5*atan2(s2,s1); % -90 to 90
% xi=0.5*atan2(s3,(sqrt((s1.^2)+ (s2.^2)))); % -45 to 45
% ximod=[-pi/2*ones(size(xi,1),1),xi,pi/2*ones(size(xi,1),1)];

%%%%  s5 display
figure(1)
%figure
subplot(1,3,1);imagesc(s0);axis image; axis off;colormap jet; title('s0');
subplot(1,3,2);imagesc(s1);axis image; axis off; title('s1');
subplot(1,3,3);imagesc(s2);axis image; axis off; title('s2');
% subplot(1,3,3);imagesc(s3); axis image;axis off;  title('s3');
saveas(gcf,'S1-3_re.png')
figure(2)
subplot(1,2,1);imagesc(shi);axis image; axis off;colormap jet; title('shi');
% subplot(1,2,2);imagesc(ximod);axis image; axis off; title('xi');
saveas(gcf,'shi and xi_re.png')
figure(3)
f_stokes_plot(s1,s2,50);



 


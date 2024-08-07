close all;clc;clear;

%%%% s1 import data to matlab from folder
I00 = importdata('I00_0001.ascii.csv');
I90 = importdata('I90_0001.ascii.csv');
I45 = importdata('I45_0001.ascii.csv'); 
I135 = importdata('I135_0001.ascii.csv');
I45q = importdata('I45q_0001.ascii.csv');
I135q = importdata('I135q_0001.ascii.csv');
imagesc(I45);axis image; axis off;colormap jet; title('s1');
%%%% s2 range of image 0 to 255
range=255;
mx=max(max(I00+I90)); 

I00=(I00/mx)*range;
I90=(I90/mx)*range;
I45=(I45/mx)*range;
I135=(I135/mx)*range;
I45q=(I45q/mx)*range;
I135q=(I135q/mx)*range;

%%%%%%%  s3 stokes parameter
s0=(I00+I90);
s1=(I00-I90)./(I00+I90);
s2=(I45-I135)./(I135+I45);
s3=(I45q-I135q)./(I135q+I45q);

%%%%%%  s4 pol angles
shi=0.5*atan2(s2,s1); % -90 to 90
xi=0.5*atan2(s3,(sqrt((s1.^2)+ (s2.^2)))); % -45 to 45
ximod=[-pi/2*ones(size(xi,1),1),xi,pi/2*ones(size(xi,1),1)];

%%%%  s5 display
figure,
subplot(1,3,1);imagesc(s1);axis image; axis off;colormap jet; title('s1');
subplot(1,3,2);imagesc(s2);axis image; axis off; title('s2');
subplot(1,3,3);imagesc(s3); axis image;axis off;  title('s3');
saveas(gcf,'S1-3_re.png')

figure,
subplot(1,2,1);imagesc(shi);axis image; axis off;colormap jet; title('shi');
subplot(1,2,2);imagesc(ximod);axis image; axis off; title('xi');
saveas(gcf,'shi and xi_re.png')

f_stokes_plot(s1,s2,s3,50);


function[Eux,Euy]=f_stokes_plot(Ue,Ve,We,n)
Ue=imresize(Ue,[n,n]);
Ve=imresize(Ve,[n,n]);
We=imresize(We,[n,n]);
Z=zeros(size(Ue));
figure; 
Q=quiver3(Z,Ue,Ve,We,2,'linewidth',1.5);
view(16,82); 
axis tight;
mags = reshape(Q.WData, numel(Q.UData), []);
currentColormap = colormap(jet);
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
set(Q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
set(Q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
 end
 


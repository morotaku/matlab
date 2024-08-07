
clear;close all; clc;

% coordinates
L=4; N = 100;
x=linspace(-L,L,N);
y=linspace(-L,L,N);
[X,Y]=meshgrid(x,-y);

% light field
lambda=530e-6;
k=2*pi/lambda;
w=1.5;
z=-1000; 
c1 = LGls(0,1,w,z,lambda,X,Y);
c2 = LGls(0,0,w,z,lambda,X,Y);
EL=0.5*(c1+c2);
ER = LGls(1,0,w,z,lambda,X,Y);

% change basis
[Ex,Ey]=f_qwp(ER,EL,45);

% stokes paramters
[s0,s1,s2,s3]=f_stokes2022(Ex,Ey);

% save picture of stokes paramters
figure(1)
subplot(1,3,1);imagesc(s1);axis image; axis off;colormap jet;
subplot(1,3,2);imagesc(s2);axis image; axis off;
subplot(1,3,3);imagesc(s3); axis image;axis off; 


% 3d plot of stokes parameters
% f_stokes_plot(X,Y,s1,s2,s3,30);

phi_12=atan2(s2,s1);
phi_23=atan2(s3,s2);
phi_31=atan2(s1,s3);
figure(2)
colormap('hsv')
subplot(2,3,1);imagesc(phi_12);axis image; axis off;colormap jet;
subplot(2,3,2);imagesc(phi_23);axis image; axis off;
subplot(2,3,3);imagesc(phi_31);axis image; axis off;

%%%%%%%%%% %%%%%%%%%% %%%%%%%%%% 
%%%%%%%%%%  funtions   %%%%%%%%%
function[s0,s1,s2,s3]=f_stokes2022(Ex,Ey)
s0=((Ex.*conj(Ex))+(Ey.*conj(Ey)));
s1=((Ex.*conj(Ex))-(Ey.*conj(Ey)));
s2=((Ex.*conj(Ey))+(Ey.*conj(Ex)));
s3=(((Ex.*conj(Ey))-(Ey.*conj(Ex))))*1i;
s1=s1./s0;s2=s2./s0;s3=s3./s0; 
s3=-s3; %%%%% doubt

end


function[Eh,Ev]=f_qwp(Ehorizontal,Evertical,sa)
QWP_h= (1)*[1,0;0,exp(1i*-pi/2)];% phase delay along horizontal
A=deg2rad(sa); % rotation of fast axis wrt to horizontal
RM = [cos(A) -sin(A); sin(A) cos(A)];
MRM = [cos(-A) -sin(-A); sin(-A) cos(-A)];
QWP=1*(MRM*QWP_h*RM );
Eh=Ehorizontal.*QWP(1,1)+Evertical.*QWP(1,2);
Ev=Ehorizontal.*QWP(2,1)+Evertical.*QWP(2,2);

end


 function[uout]=LGls(l,p,w0,z,lambda,X,Y)

phi=angle(X+1i*Y);rho=sqrt(X.^2+Y.^2);
k=2*pi/lambda;zr=pi*w0^2/lambda;
w=w0*sqrt(1+(z/zr)^2);
R=z*(1+(zr/z)^2);
La=Laguerre(p,abs(l),2*rho.^2/(w^2));
LG=w0/w*sqrt(2*factorial(p)/(pi*factorial(abs(l)+p)))*(2*rho.^2/w^2).^(abs(l)/2).*La.*exp(1i*(2*p+l+1)*atan(z/zr)).*exp(-rho.^2/w^2).*exp(-1i*k*rho.^2/R).*exp(1i*l*phi);    
uout=LG/sqrt(sum(sum(abs(LG).^2)));
  end
 

 function[out]=circ(r)
out=abs(r)<=1;
 end


 function y=Laguerre(p,l,x)
y=zeros(p+1,1);
if p==0
    y=1;
else
for m=0:p
    y(p+1-m)=((-1).^m.*(factorial(p+l)))./(factorial(p-m).*factorial(l+m).*factorial(m));
end
end
y=polyval(y,x);
 end
 
  function[Eux,Euy]=f_stokes_plot(X,Y,Ue,Ve,We,n)

X=imresize(X,[n,n]);
Y=imresize(Y,[n,n]);
Ue=imresize(Ue,[n,n]);
Ve=imresize(Ve,[n,n]);
We=imresize(We,[n,n]);

Z=zeros(size(X));
figure; 
Q=quiver3(X,Y,Z,Ue,Ve,We,2,'linewidth',1.5);
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
saveas(gcf,'z stokes','png')
% close all;
 end
 
 



 
 
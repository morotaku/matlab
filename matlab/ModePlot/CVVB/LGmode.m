function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l).*laguerre(p, l, 2*r.^2/W^2).*exp(Gouy);
LG=abs(LG1.^2);
y=LG./max(LG(:));
end

function y=laguerre(n,l,x)

%LAGUERRE Laguerre polynomial.
if numel(n)>1
    error('Laser Toolbox: HERMITE: N should be a non-negative scalar integer.');
end

if (n<0 || abs(round(n)-n)>0)
    abs(round(n)-n)>0
   error('Laser Toolbox: HERMITE: N should be a non-negative integer.');
end

if prod(size(l))>1
    error('Laser Toolbox: HERMITE: L should be a non-negative scalar integer.');
end

if (l<0 || abs(round(l)-l)>0)
    abs(round(n)-n)>0
   error('Laser Toolbox: HERMITE: L should be a non-negative integer.');
end

if n==0
   y=1;
elseif n==1
   y=-1*x+1+l;
elseif n>1
   y= (2*n+l-1-x)./n.*laguerre(n-1,l,x) - (n+l-1)/n*laguerre(n-2,l,x);
end
end

function y=hgb(p,l,r,z,w,k,A,B,C,D,G0)
i=1j;
C1=(i*k*G0*factorial(p)/(2*B*w^(2*p)))*((1/(w^2)+(i*k*A/(2*B)))^(-1-p))*exp(-i*k*z);
C2=exp(((-i*k*D.*(r.^2))./(2*B)));
C3=exp(-((k.*r)./(2*B)).^2./(1/w^2+(i*k*A/(2*B))));
lag=laguerre(p,l,(k/(2*B))^2/((1/w^2)+(i*k*A/(2*B))).*(r.^2));
hgb=C1.*C2.*C3.*lag;
HGM=abs(hgb.^2);
y=HGM./max(HGM(:));
end
function y=f_hgb(n,s,f,z,lam,w,w2,G0,r)
    %ABCD matrix
    Z=z+s+f;
    A=1+(s-Z)/f;
    B=Z+(s^2-2*s)/f;
    C=-1/f;
    D=1-s/f;
    k=2*pi/lam; %Wave number
    zr=w^2*k/2; %Rayleigh length
    q0=1j*k*w^2/2;
    Q=(A+B/q0)/(C+D/q0);
    W=w*sqrt(A^2+(B*lam./(pi*w^2))^2);

    C1=G0*factorial(n)/2^n;
    C2=0;
    %A=(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1))*exp(-1j*k*z).*laguerre(m,2/w^2.*r.^2)%.*exp(-1j*k/(2*Q).*r.^2)
    for m=0:n
        C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1))*exp(-1j*k*Z).*laguerre(m,2/w2^2.*r.^2).*exp(-1j*k/(2*Q).*r.^2);
        %C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1)).*laguerre(m,2/w^2.*r.^2).*exp(-1j*k/(2*Q).*r.^2);
    end
    
    lag=C1.*C2;
    HGM=abs(lag).^2;
    y=HGM/max(HGM(:));
end
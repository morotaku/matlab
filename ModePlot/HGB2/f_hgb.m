function y=f_hgb(n,f,z,lam,w0,G0,r)
    %ABCD matrix
    s=0;
    A=-z./f;
    B=-(z./f).*s+f+z;
    C=-1/f;
    D=1-s/f;
    k=2*pi/lam; %Wave number
    q0=1j*k.*w0.^2./2;
    %Q=(A+B./q0)./(C+D./q0);
    %W=w0.*sqrt(A.^2+(B.*lam./(pi*w0^2)).^2);

    %%%HGB
    C1=1j*k*G0*factorial(n)./(2*w0^(2*n).*B).*((1./w0.^2)+(1j*k.*A./(2.*B))).^(-1-n).*exp(-1j*k.*z);
    C2=exp(-1j*k.*(D.*r.^2)./(2.*B));
    C3=exp(-(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B)))).*laguerre(n,(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B))));
    hgb=C1.*C2.*C3;
    hgb_n=abs(hgb)./max(max(abs(hgb)));
    HGB=hgb_n.^2;
    HGB_n=HGB./max(HGB(:));
    y=HGB_n;
end

function y=laguerre(n,x)
y=0;
    if n==0
        y=1;
    else
        for m=0:n
            y=y+(factorial(n)^2*(-1)^m.*x.^m)./(factorial(m)*factorial(m)*factorial(n-m));
        end
    end
end
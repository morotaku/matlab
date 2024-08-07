function y=HGmode_45(n,x,y,w,r,z,R,k)
    hg=0;
    for i=0:n
        if rem(i,2)==1
            hg=hg-nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        else
            hg=hg+nchoosek(n,i).*HGmode(i,n-i,x,y,w,r,z,R,k);
        end
    end
    y=hg./(sqrt(2)^n);
    %y=hg./max(max(abs(hg)));
end

function y=HGmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-(r.^2./W^2)-(1j*k.*r.^2/(2*R)));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;
    %y=u./max(max(abs(u)));
end


function y=HermitePol(n,x)
Hn=0;
for m=0:floor(n/2)
    Hn=Hn+((-1)^m.*(2.*x).^(n-2*m))./(factorial(m)*factorial(n-2*m));
end
Hn=Hn.*factorial(n);
y=Hn;
end
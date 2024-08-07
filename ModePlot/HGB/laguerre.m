function y=laguerre(n,x)
    if n==0
        y=1;
    else
        y=zeros(size(x));
        for k=0:n
            y=y+(-1)^k*factorial(n).*x.^k./(factorial(n-k)*(factorial(k))^2);
        end
        y=factorial(n).*y;
    end
end
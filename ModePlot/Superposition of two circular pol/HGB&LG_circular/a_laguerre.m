function y=a_laguerre(n,x)
y=0;
    if n==0
        y=1;
    elseif n==1
        
        for m=0:n
            y=y+(factorial(n)^2*(-1)^m.*x.^m)./(factorial(m)*factorial(m+k)*factorial(n-m-k));
        end
    end
end
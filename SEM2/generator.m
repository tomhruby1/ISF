function y = generator(a,b,c,u,lambda,N)
   
    y = zeros(1,N);
    na = length(a);
    nc = length(c);
    e = lambda*randn(1,N+nc-1);

    y(1) = e(1);
    for i = 2:N
        y(i) = -a(1) * y(i-1) + b(1)*u(i-1) + e(i) + c*e(i-1);
    end  
end
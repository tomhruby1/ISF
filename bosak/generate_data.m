function y = generate_data(a, b, c, u, lambda, N);
    na = length(a);
    nb = length(b);
    nc = length(c);

    e = randn(N+nc-1, 1)*lambda;
    
    y = zeros(N+1, 1);
%     bacha t = t-1
    for i = 2:1:N+1
        y(i, 1) = c*e(i+nc-2:-1:i-1, 1);
        if na > 1
            if i > 1 && i-na+1 < 1 %mame jen nektera data z minulosti                                     
                sumY = -a(2:1:i)*y(i-1:-1:1, 1);
            elseif  i-1 >= 1 && i-na+1 >= 1 %mame vsechna potrebna data
                sumY = -a(2:1:na)*y(i-1:-1:i-na+1, 1);          
            end
        end
        sumU = 0;
            if i-nb < 1 %mame jen nektera data z minulosti                                     
                sumU = b(1:i-1)*u(i-1:-1:1, 1);
            elseif  i-1 >= 1 && i-nb+1 >= 1 %mame vsechna potrebna data
                sumU = b(1:1:nb)*u(i-1:-1:i-nb, 1);          
            end
        
        y(i, 1) = y(i, 1) + sumY + sumU;
                                    
    end             
end
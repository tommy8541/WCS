B = [0.01 0.03 0.05 0.1]; %blocking rate
m = [120 60 40];   %channel number
rho = zeros(length(m), length(B));
for i = 1:length(m)
    for j = 1:length(B)
        rho(i,j) = rho_operator(m(i),B(j));
    end
end

function mean = rho_operator(m, B)  %use dichotomy
format long
%B = 0.01;
%m = 120;
LS = 0;
RS = 2*m;
E = 1;
syms k;


while(E > 0.0000001)
    mean = (LS + RS)/2;
    %block_rate = mean^m/(expand(factorial(sym(m)))*symsum(mean^k/expand(factorial(sym(k))),k,0,m));
    sum = 0;
    for k = 0:m
       sum = sum + mean^(k)/prod(vpa(1:k,50));
    end
    block_rate = mean^(m)/(sum*prod(vpa(1:m,50))); 
    %block_rate = mean^m/(prod(vpa(1:m, 50))*symsum(mean^k/prod(vpa(1:k, 50)),k,0,m));
    E = abs(block_rate - B);
    if(block_rate > B)
        RS = mean;
    else
        LS = mean;
    end
end
end



%disp(mean);
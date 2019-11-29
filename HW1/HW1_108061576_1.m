rho = zeros(41, 4);
B = [0.01 0.03 0.05 0.1];
m = [1:20 200:220];
for i = 1:length(m)
    for j= 1:length(B)
        rho(i, j) = find_rho (m(i), B(j));
    end
end

function mean = find_rho(m, B)
format long
%B = 0.01;
%m = 20;
LS = 0;
RS = 2*m;
syms k;

count = 0;

while(count < 50)
    count = count + 1;
    %disp(count);
    mean = (LS + RS)/2;
    block_rate = sym(mean)^m / (factorial(sym(m))*symsum(sym(mean)^k/factorial(sym(k)),k,0,m));
    %sum = 0;
    %for k = 0:m
    %   sum = sum + sym(mean)^(k)/factorial(sym(k));
    %end
    %block_rate = sym(mean)^(m)/(sum*factorial(sym(m))); 
    %E = abs(log(block_rate) - log(B));
    if(log(block_rate) > log(B))
        RS = mean;
    else
        LS = mean;
    end
end
end


%disp(mean);
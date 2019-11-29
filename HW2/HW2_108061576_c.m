%% Initialization
clear all;
close all;
clc;

v = unifrnd(20, 90, [1 2*10^6]); %km/hr
theta = unifrnd(-pi, pi, [1 2*10^6]);  %uniformly distributed angle
v_change = v .* (1000/3600);
fc = 2 * (10 ^ 9 ); %Hz
fm = v_change * fc / ( 3 * (10^8));

fD = cos(theta) .* fm;

fD_min = min(fD);
fD_max = max(fD);
[countsfD, binsfD] = hist(fD, 10^3);
cdf = cumsum(countsfD) / sum(countsfD);

pdf = zeros(1,length(cdf));
for i = 1:length(cdf)
    if i ==1
        pdf(i) = cdf(i) / ((fD_max-fD_min)/(1000-1));
    else
        pdf(i) = (cdf(i) - cdf(i-1)) / ((fD_max-fD_min)/(1000-1));
    end
end
%pdf = diff(cdf);
xs = fD_min:(fD_max-fD_min)/((10^3)-1):fD_max;
figure(1)
plot(xs, cdf); 
title('Cumulative Distribution Function v.s. Doppler Shift');
xlabel('Doppler Shift');
ylabel('CDF');
print('-f1', '-djpeg', '-r300', 'hw2_c_cdf');

figure(2)
%pdf =smooth(diff(cdf));
plot(xs, pdf);
title("Probability Density Funcion v.s. Doppler Shift");
xlabel("Doppler Shift");
ylabel("PDF");
print('-f2', '-djpeg', '-r300', 'hw2_c_pdf');



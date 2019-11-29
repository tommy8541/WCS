%% Initialization
clear all;
close all;
clc;

v_a = 20; %km/hr
v_a_change = v_a*1000/3600;
fc_a = 2 * (10 ^ 9 ); %Hz
fm_a = v_a_change * fc_a / ( 3 * (10^8));

v_b = 90; %km/hr
v_b_change = v_b*1000/3600;
fc_b = 26 * (10 ^ 9 ); %Hz
fm_b = v_b_change * fc_b / ( 3 * (10^8));

gamma = -1:0.001:1;
gamma_pdf_a = gamma * fm_a;
gamma_pdf_b = gamma * fm_b;

%pdf
pdf_a = 1 ./ ( pi .*sqrt(1- gamma.^2)) ;

%cdf
%cdf = acos(-gamma)/pi ;
cdf_a = 1 - acos(gamma) ./ pi;

pdf_b = 1 ./ ( pi .*sqrt(1- gamma.^2)) ;
cdf_b = 1 - acos(gamma) ./ pi;

figure(1)
plot(gamma_pdf_a, pdf_a)
title("Theoretical results of PDF for a");
xlabel("Normalized Doppler Shift");
ylabel("PDF of a");
print('-f1', '-djpeg', '-r300', 'hw2_d_pdf_a');

figure(2)
plot(gamma_pdf_a, cdf_a)
title("Theoretical results of CDF for a");
xlabel("Normalized Doppler Shift");
ylabel("CDF of a");
print('-f2', '-djpeg', '-r300', 'hw2_d_cdf_a');

figure(3)
plot(gamma_pdf_b, pdf_b)
title("Theoretical results of PDF for b");
xlabel("Normalized Doppler Shift");
ylabel("PDF_b");
print('-f3', '-djpeg', '-r300', 'hw2_d_pdf_b');

figure(4)
plot(gamma_pdf_b, cdf_b)
title("Theoretical results of CDF for b");
xlabel("Normalized Doppler Shift");
ylabel("CDF_b");
print('-f4', '-djpeg', '-r300', 'hw2_d_cdf_b');

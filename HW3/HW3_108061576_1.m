clear all;
close all;

fmT = [0.01, 0.1, 0.5];
length = 10000;
t = 1:300;
gi = zeros(3, length);
gq = zeros(3, length);
zeta = 2 - cos(pi .* fmT ./ 2) - sqrt(( 2 - cos(pi .* fmT ./ 2)).^2 - 1);

sigma_001_1 = normrnd(0, sqrt((1+zeta(1))/(1-zeta(1)/2)), 1 ,length);
sigma_01_1 = normrnd(0, sqrt((1+zeta(2))/(1-zeta(2)/2)), 1 ,length);
sigma_05_1 = normrnd(0, sqrt((1+zeta(3))/(1-zeta(3)/2)), 1 ,length);

sigma_001_2 = normrnd(0, sqrt((1+zeta(1))/(1-zeta(1)/2)), 1 ,length);
sigma_01_2 = normrnd(0, sqrt((1+zeta(2))/(1-zeta(2)/2)), 1 ,length);
sigma_05_2 = normrnd(0, sqrt((1+zeta(3))/(1-zeta(3)/2)), 1 ,length);

%IV for gi
gi(1,1) = (1-zeta(1)) * sigma_001_1(1);
gi(2,1) = (1-zeta(2)) * sigma_01_1(1);
gi(3,1) = (1-zeta(3)) * sigma_05_1(1);

%IV for gq
gq(1,1) = (1-zeta(1)) * sigma_001_2(1);
gq(2,1) = (1-zeta(2)) * sigma_01_2(1);
gq(3,1) = (1-zeta(3)) * sigma_05_2(1);

for i = 2:length
    gi(1, i) = zeta(1) .* gi(1, i-1) + (1-zeta(1)) .* sigma_001_1(1, i-1);
    gi(2, i) = zeta(2) .* gi(2, i-1) + (1-zeta(2)) .* sigma_01_1(1, i-1);
    gi(3, i) = zeta(3) .* gi(3, i-1) + (1-zeta(3)) .* sigma_05_1(1, i-1);

    gq(1, i) = zeta(1) .* gq(1, i-1) + (1-zeta(1)) .*sigma_001_2(1, i-1);
    gq(2, i) = zeta(2) .* gq(2, i-1) + (1-zeta(2)) .*sigma_01_2(1, i-1);
    gq(3, i) = zeta(3) .* gq(3, i-1) + (1-zeta(3)) .*sigma_05_2(1, i-1);
end
 
r = zeros(3, length);
for i = 1:length
    r(1, i) = gi(1, i) + 1j * gq(1, i);
    r(2, i) = gi(2, i) + 1j * gq(2, i);
    r(3, i) = gi(3, i) + 1j * gq(3, i);
end

r_head = zeros(3, length);
for i = 1:length
    r_head(1, i) = 5*log10(gi(1, i)^2 + gq(1, i)^2);
    r_head(2, i) = 5*log10(gi(2, i)^2 + gq(2, i)^2);
    r_head(3, i) = 5*log10(gi(3, i)^2 + gq(3, i)^2);
end

figure(1)
subplot(3, 1, 1);
plot(t, r_head(1, 1:300));
title("channel output for fmT = 0.01")
xlabel("t/T")
ylabel("Envelope Level (dB)")

subplot(3, 1, 2);
plot(t, r_head(2, 1:300));
title("channel output for fmT = 0.1")
xlabel("t/T")
ylabel("Envelope Level (dB)")

subplot(3, 1, 3);
plot(t, r_head(3, 1:300));
title("channel output for fmT = 0.5")
xlabel("t/T")
ylabel("Envelope Level (dB)")
print('-f1', '-djpeg', '-r300', 'channel output');   

%autocorrelation
acf_1 = autocorr(r(1 , :), 10/fmT(1));
acf_2 = autocorr(r(2 , :), 10/fmT(2));
acf_3 = autocorr(r(3 , :), 10/fmT(3));

figure(2)
plot(0:fmT(1):10, acf_1);
hold on
plot(0:fmT(2):10, acf_2, '-r');
hold on
plot(0:fmT(3):10, acf_3, '-k');
legend('fmT = 0.01','fmT = 0.1','fmT = 0.5');
title("channel output autocorrelation ")
xlabel("fm£n")
ylabel("Autocorrelation")
print('-f2', '-djpeg', '-r300', 'autocorrelation');   



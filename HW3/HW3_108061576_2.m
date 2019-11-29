clear all;
close all;

M = [8, 16];
N = 4 * M + 2;
fmT = [0.01, 0.1, 0.5];
length = 10000;
alpha = 0;
t = 1:300;

beta = zeros(2, 16);
for i = 1:M(1)
    beta(1,i) = (pi / M(1)) * i;
end

for i = 1:M(2)
    beta(2,i) = (pi / M(2)) * i;
end

sum_I_8 = zeros(3, length);
sum_Q_8 = zeros(3, length);

% M = 8
for j = 1:10000
    sum_I_8(1,length) = 0;
    sum_Q_8(1,length) = 0;
    for i = 1:M(1)
        sum_I_8(1,j) = cos(beta(1,i))*cos(2*pi*fmT(1)*cos(2*pi*i/N(1))*j) + sum_I_8(1,j);
        sum_Q_8(1,j) = sin(beta(1,i))*cos(2*pi*fmT(1)*cos(2*pi*i/N(1))*j) + sum_Q_8(1,j);
    end
end

for j = 1:10000
    sum_I_8(1,length) = 0;
    sum_Q_8(1,length) = 0;    
    for i = 1:M(1)
        sum_I_8(2,j) = cos(beta(1,i))*cos(2*pi*fmT(2)*cos(2*pi*i/N(1))*j) + sum_I_8(2,j);
        sum_Q_8(2,j) = sin(beta(1,i))*cos(2*pi*fmT(2)*cos(2*pi*i/N(1))*j) + sum_Q_8(2,j);
    end
end

for j = 1:10000
    sum_I_8(1,length) = 0;
    sum_Q_8(1,length) = 0;    
    for i = 1:M(1)
        sum_I_8(3,j) = cos(beta(1,i))*cos(2*pi*fmT(3)*cos(2*pi*i/N(1))*j) + sum_I_8(3,j);
        sum_Q_8(3,j) = sin(beta(1,i))*cos(2*pi*fmT(3)*cos(2*pi*i/N(1))*j) + sum_Q_8(3,j);
    end
end

gi_8 = zeros(3, length);
gq_8 = zeros(3, length);
for i = 1:3
    for j = 1:10000
        gi_8(i,j) = sqrt(2) * (2 * sum_I_8(i,j) + sqrt(2)*cos(alpha)*cos(2*pi*fmT(i)*j));
        gq_8(i,j) = sqrt(2) * (2 * sum_Q_8(i,j) + sqrt(2)*sin(alpha)*cos(2*pi*fmT(i)*j));
    end
end

g_8 = zeros(3, length);
for i = 1:10000
    g_8(1,i) = gi_8(1, i) + 1j * gq_8(1,i);
    g_8(2,i) = gi_8(2, i) + 1j * gq_8(2,i);
    g_8(3,i) = gi_8(3, i) + 1j * gq_8(3,i);
end

g_head_8 = zeros(3, length);
for i = 1:10000
    g_head_8(1,i) = 5*log10(gi_8(1, i)^2 + gq_8(1, i)^2);
    g_head_8(2,i) = 5*log10(gi_8(2, i)^2 + gq_8(2, i)^2);
    g_head_8(3,i) = 5*log10(gi_8(3, i)^2 + gq_8(3, i)^2);
end

figure(1)
subplot(3, 1, 1);
plot(t, g_head_8(1, 1:300));
title("M = 8 channel output for fmT = 0.01")
xlabel("t/T")
ylabel("Envelope Level (dB)")

subplot(3, 1, 2);
plot(t, g_head_8(2, 1:300));
title("M = 8 channel output for fmT = 0.1")
xlabel("t/T")
ylabel("Envelope Level (dB)")

subplot(3, 1, 3);
plot(t, g_head_8(3, 1:300));
title("M = 8 channel output for fmT = 0.5")
xlabel("t/T")
ylabel("Envelope Level (dB)")
print('-f1', '-djpeg', '-r300', 'M=8 channel output');   
    
% M = 16
sum_I_16 = zeros(3, length);
sum_Q_16 = zeros(3, length);
for j = 1:10000
    sum_I_16(1,length) = 0;
    sum_Q_16(1,length) = 0;
    for i = 1:M(2)
        sum_I_16(1,j) = cos(beta(2,i))*cos(2*pi*fmT(1)*cos(2*pi*i/N(2))*j) + sum_I_16(1,j);
        sum_Q_16(1,j) = sin(beta(2,i))*cos(2*pi*fmT(1)*cos(2*pi*i/N(2))*j) + sum_Q_16(1,j);
    end
end

for j = 1:10000
    sum_I_16(1,length) = 0;
    sum_Q_16(1,length) = 0;    
    for i = 1:M(2)
        sum_I_16(2,j) = cos(beta(2,i))*cos(2*pi*fmT(2)*cos(2*pi*i/N(2))*j) + sum_I_16(2,j);
        sum_Q_16(2,j) = sin(beta(2,i))*cos(2*pi*fmT(2)*cos(2*pi*i/N(2))*j) + sum_Q_16(2,j);
    end
end

for j = 1:10000
    sum_I_16(1,length) = 0;
    sum_Q_16(1,length) = 0;    
    for i = 1:M(2)
        sum_I_16(3,j) = cos(beta(2,i))*cos(2*pi*fmT(3)*cos(2*pi*i/N(2))*j) + sum_I_16(3,j);
        sum_Q_16(3,j) = sin(beta(2,i))*cos(2*pi*fmT(3)*cos(2*pi*i/N(2))*j) + sum_Q_16(3,j);
    end
end

gi_16 = zeros(3, length);
gq_16 = zeros(3, length);
for i = 1:3
    for j = 1:10000
        gi_16(i,j) = sqrt(2) * (2 * sum_I_16(i,j) + sqrt(2)*cos(alpha)*cos(2*pi*fmT(i)*j));
        gq_16(i,j) = sqrt(2) * (2 * sum_Q_16(i,j) + sqrt(2)*sin(alpha)*cos(2*pi*fmT(i)*j));
    end
end

g_16 = zeros(3, length);
for i = 1:10000
    g_16(1,i) = gi_16(1, i) + 1j * gq_16(1,i);
    g_16(2,i) = gi_16(2, i) + 1j * gq_16(2,i);
    g_16(3,i) = gi_16(3, i) + 1j * gq_16(3,i);
end

g_head_16 = zeros(3, length);
for i = 1:10000
    g_head_16(1,i) = 5*log10(gi_16(1, i)^2 + gq_16(1, i)^2);
    g_head_16(2,i) = 5*log10(gi_16(2, i)^2 + gq_16(2, i)^2);
    g_head_16(3,i) = 5*log10(gi_16(3, i)^2 + gq_16(3, i)^2);
end

figure(2)
subplot(3, 1, 1);
plot(t, g_head_16(1, 1:300));
title("M = 16 channel output for fmT = 0.01")
xlabel("t/T")
ylabel("Envelope Level (dB)")

subplot(3, 1, 2);
plot(t, g_head_16(2, 1:300));
title("M = 16 channel output for fmT = 0.1")
xlabel("t/T")
ylabel("Envelope Level (dB)")

subplot(3, 1, 3);
plot(t, g_head_16(3, 1:300));
title("M = 16 channel output for fmT = 0.5")
xlabel("t/T")
ylabel("Envelope Level (dB)")
print('-f2', '-djpeg', '-r300', 'M=16 channel output');   

acf_8_1 = autocorr(g_8(1,:), 10/fmT(1));
acf_8_2 = autocorr(g_8(2,:), 10/fmT(2));
acf_8_3 = autocorr(g_8(3,:), 10/fmT(3));

figure(3)
plot(0:fmT(1):10, acf_8_1);
hold on
plot(0:fmT(2):10, acf_8_2, '-r');
hold on
plot(0:fmT(3):10, acf_8_3, '-k');
legend('fmT = 0.01','fmT = 0.1','fmT = 0.5');
title("M = 8, channel output autocorrelation")
xlabel("fm£n")
ylabel("Autocorrelation")
print('-f3', '-djpeg', '-r300', 'M=8 autocorr');

acf_16_1 = autocorr(g_16(1,:), 10/fmT(1));
acf_16_2 = autocorr(g_16(2,:), 10/fmT(2));
acf_16_3 = autocorr(g_16(3,:), 10/fmT(3));

figure(4)
plot(0:fmT(1):10, acf_16_1);
hold on
plot(0:fmT(2):10, acf_16_2, '-r');
hold on
plot(0:fmT(3):10, acf_16_3, '-k');
legend('fmT = 0.01','fmT = 0.1','fmT = 0.5');
title("M = 16, channel output autocorrelation")
xlabel("fm£n")
ylabel("Autocorrelation")
print('-f4', '-djpeg', '-r300', 'M=16 autocorr');
f   = 40;
fs  = 1000;
c   = 10;

pm  = morlet_test(f, c, fs);
nm  = morlet_test(-f, c, fs);
spm = fft(pm);
snm = fft(nm);
ppm = abs(spm) .^ 2;
pnm = abs(snm) .^ 2;

t   = floor(numel(pm)) / fs;
t   = linspace(-t, t, numel(pm));
f   = linspace(-500, 500, numel(spm));

figure;

subplot(3, 2, 1);
title('FPos');
xlabel('Time (s)');
ylabel('Amplitude');
hold on;
plot(t, real(pm));
plot(t, imag(pm));

subplot(3, 2, 2);
title('FNeg');
xlabel('Time (s)');
ylabel('Amplitude');
hold on;
plot(t, real(nm));
plot(t, imag(nm));

subplot(3, 2, 3);
title('S(FPos)');
xlabel('Frequency (Hz)');
ylabel('Fourier coefficient');
hold on;
plot(f, real(spm));
plot(f, imag(spm));

subplot(3, 2, 4);
title('S(FNeg)');
xlabel('Frequency (Hz)');
ylabel('Fourier coefficient');
hold on;
plot(f, real(snm));
plot(f, imag(snm));

subplot(3, 2, 5);
title('P(FPos)');
xlabel('Frequency (Hz)');
ylabel('Fourier coefficient');
hold on;
plot(f, ppm);

subplot(3, 2, 6);
title('P(FNeg)');
xlabel('Frequency (Hz)');
ylabel('Fourier coefficient');
hold on;
plot(f, pnm);



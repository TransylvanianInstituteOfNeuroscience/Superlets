clear all;
close all;

fs  = 1024;
t   = 0:(1/fs):8;
s   = chirp(t, 20, t(end), 380);

% compute spectrogram and scalograms
f   = 1:0.5:400;
o   = [10, 40];
mul = 1; % enable multiplicative super-resolution (for better frequency localization)

asc = aslt(s, fs, f, 3, o, mul);
fsc = faslt(s, fs, f, 3, o, mul);
[ft, ft_f, ft_t] = spectrogram(s, gausswin(102), 101, 2048, fs);
ft = abs(ft).^2; % compute the real-valued power spectrum


figure;
subplot(1,3,1);
imagesc(t, f, asc);
set(gca, 'ydir', 'normal');
title('ASLT');
subplot(1,3,2);
imagesc(t, f, fsc);
set(gca, 'ydir', 'normal');
title('FASLT');
subplot(1,3,3);
imagesc(ft_t, ft_f, ft);
set(gca, 'ydir', 'normal');
title('FOURIER');
clear all
close all
load('rngpro.mat');

fs  = 1024;
s   = sx_1;%(fs * 20 : fs * 30);
f   = [-100, 100];
nf  = 101; % 2 hz per bin

% compute scalogram
axt = (1 / fs) : (1 / fs) : (numel(s) / fs);
axf = linspace(f(1), f(2), nf);
sl  = nfaslt(s, fs, f, nf, 3, [1, 20], 1);

% I found the log gives a cleaner-looking result
sl = log10(sl + 0.0001);

% plot
figure;
imagesc(axt, axf, sl);
set(gca, 'ydir', 'normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
colormap jet;
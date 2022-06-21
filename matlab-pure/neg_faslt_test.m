
% define frequencies
fs      = 1000;
f       = 10;
f_pos   = 10;
f_neg   = -10;
range   = [-20, 20];
n_bins  = 41;

% create signals
t       = 0:(1/fs):1;
s       = sin(2 * pi * f * t);
s_pos   = exp(2i * pi * f_pos * t);
s_neg   = exp(2i * pi * f_neg * t);

% compute superlet
spec        = nfaslt(s,     fs, [-20, 20], 41, 3, [1, 10], 0);
spec_pos    = nfaslt(s_pos, fs, [-20, 20], 41, 3, [1, 10], 0); 
spec_neg    = nfaslt(s_neg, fs, [-20, 20], 41, 3, [1, 10], 0); 

% plot
subplot(1, 3, 1);
hold on;
imagesc(t, linspace(range(1), range(2), n_bins), spec);
set(gca, 'ydir', 'normal');
title('Real plane wave');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(1, 3, 2);
hold on;
imagesc(t, linspace(range(1), range(2), n_bins), spec_pos);
set(gca, 'ydir', 'normal');
title('Positive-frequency plane wave');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(1, 3, 3);
hold on;
imagesc(t, linspace(range(1), range(2), n_bins), spec_neg);
set(gca, 'ydir', 'normal');
title('Negative-frequency plane wave');
xlabel('Time (s)');
ylabel('Frequency (Hz)');









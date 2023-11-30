
% params
n   = 1000;
fs  = 1000;
fr  = [10 90];

c1  = 5;
sr  = [5 20];

% create signal
t   = (1 : n) / fs;

% s   = sin(2 * pi * 40 * t);               % one sine
s   = 3 * sin(2 * pi * 40 * t) + ...
      3 * sin(2 * pi * 70 * t);           % two sines (single trial)
% s   = [3 * sin(2 * pi * 40 * t); ...
%        3 * sin(2 * pi * 70 * t)];         % two sines (2 trials)
% s   = rand(1, n);                         % random signal
  
% compute
x   = faslt(s, fs, fr, nf, c1, sr, 1, 1);

% plot
figure;
x_axis = linspace(1/fs, n/fs, n);
y_axis = linspace(fr(1), fr(2), nf);
imagesc(x_axis, y_axis, x);
set(gca, 'ydir', 'normal');
colormap jet;
colorbar;




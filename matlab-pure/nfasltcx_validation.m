%% Detailed test script for nfasltcx.m function
%
% Validates that nfasltcx produces complex-valued output whose magnitude
% squared equals the real-valued output from nfaslt. Uses the same
% synthetic signal as the original "Superlets_toy_Data" demo script.
%
% Ben Jancovich, 2025
% Centre for Marine Science and Innovation
% School of Biological, Earth and Environmental Sciences
% University of New South Wales, Sydney, Australia
%
clear
close all
clc

%% Global and SLT parameters

fs = 1024; % Sampling frequency

fMin = 10; % Low frequency limit (Hz)
fMax = 80; % High frequency limit (Hz)
fStep = 0.5; % Frequency Step (Hz)

c1      = 3; % Number of cycles in initial superlet.
o   = [5, 10]; % Interval of superresolution orders.
mult    = 1; % Multiplicative or Additive Superresolution

% Set up frequency inputs
frequency_vector = fMin:fStep:fMax;
Fi      = [min(frequency_vector), max(frequency_vector)];
Nf      = length(frequency_vector);

% Test pass/fail criteria - Max delta between power SLTs must be < tol
tolerance = 1e-12;

%% Generate Test signal

N = fs * 3.5; % packet length
vfTarget =                  [20, 40, 60];   %target frequencies
vfNeighbF=                  [30, 50, 70];   %neighboring frequencies
nFreqs =                    numel(vfTarget);
nWaveCycles =               11;
nTNeighbSpacingInBursts=    1 + 2/nWaveCycles;
bTNeighbRelativeSpacing =   true;
nPacketSpacinginBursts =    1/nWaveCycles;
nPackets =                  2;
nPreSpaceS =                0.25;

nLongestBurstlen =  round(nWaveCycles * fs / min([vfTarget vfNeighbF]));
nPacketSpacing =    round(nPacketSpacinginBursts * nLongestBurstlen);
nPacketLen =        round(nLongestBurstlen * (2 + nTNeighbSpacingInBursts - 1));

xTrg = zeros(1,N);  %a place to store the target
xNF = zeros(1,N);   %a place to store the neighbor in frequency
xNT = zeros(1,N);   %a place to store the neighbor in time

for i = 1 : numel(vfTarget)         %for all target frequencies
    nPacketOffset = nPreSpaceS * fs + (i - 1)*(nPacketLen + nPacketSpacing);    %initial offset - useful to place at different time moments packets  
    fTarg =  vfTarget(i);
    fNeigF = vfNeighbF(i);
    nTarg =  round(nWaveCycles * (fs / fTarg));   %length of target samples
    vTarg =  sin(2*pi*fTarg/fs  * (0 : nTarg-1)); %target burst

    vfNeighF = linspace(fNeigF, fTarg, nPackets);   %neighboring frequencies
    if bTNeighbRelativeSpacing %same spacing no matter the frequency
        vfNeighT = linspace(round(nTNeighbSpacingInBursts * nTarg),0, nPackets); %dictated by the lowest frequency
    end

    for p = 1 : nPackets - 1
        fNeigF =    vfNeighF(p);
        nNeigF =    nTarg;    %length of the freq neighbor burst   
        vNeighbF =  sin(2*pi*fNeigF/fs * (0 : nNeigF-1) - pi/1.5);

        nFirstTargSample =  nPacketOffset +  (p - 1) * nFreqs * (nPacketLen + nPacketSpacing) + round((nLongestBurstlen-nTarg ) / 2);
        nFirstNFSample =    nFirstTargSample;
        nFirstNTSample =    nFirstTargSample + round(vfNeighT(p));
        
        xTrg(nFirstTargSample : nFirstTargSample + nTarg -  1) = vTarg; 
        xNF (nFirstNFSample   : nFirstNFSample +   nNeigF - 1) = vNeighbF;
        xNT (nFirstNTSample   : nFirstNTSample  +  nTarg -  1) = vTarg;
    end
end

xSignal = xTrg + xNF + xNT;

%% Run Real and Complex Superlet Transform Functions

tic;
slt_real = nfaslt(xSignal, fs, Fi, Nf, c1, o, mult);
t_real = toc;

tic;
slt_complex = nfasltcx(xSignal, fs, Fi, Nf, c1, o, mult);
t_complex = toc;

%% Extract power from complex result

slt_power_from_complex = abs(slt_complex).^2;

%% Compare results

diff_matrix = slt_real - slt_power_from_complex;
max_abs_diff = max(abs(diff_matrix(:)));
max_rel_diff = max_abs_diff / max(slt_real(:));

%% Detailed analysis

fprintf('Real function time:     %.6f seconds\n', t_real);
fprintf('Complex function time:  %.6f seconds\n', t_complex);
fprintf('Signal length: %d samples\n', length(xSignal));
fprintf('Frequency range: %.1f - %.1f Hz (%d points)\n', Fi(1), Fi(2), Nf);

% Statistics of outputs
fprintf('\nSLT statistics:\n');
fprintf('NNFASLT power = nfaslt_slt\n')
fprintf('\tMin: %.6f, Max: %.6f, Mean: %.6f\n', ...
    min(slt_real(:)), max(slt_real(:)), mean(slt_real(:)));
fprintf('NFASLTCX power = |nfasltcx_slt|^2\n')
fprintf('\tMin: %.6f, Max: %.6f, Mean: %.6f\n', ...
    min(slt_power_from_complex(:)), max(slt_power_from_complex(:)), mean(slt_power_from_complex(:)));

% Difference statistics
fprintf('\nSLT Difference statistics:\n');
fprintf('Relative differences = nfaslt_slt - |nfasltcx_slt|^2\n ')
fprintf('\tMin: %.2e, Max: %.2e, Mean: %.2e\n', ...
    min(abs(diff_matrix(:))./slt_real(:)), max(abs(diff_matrix(:))./slt_real(:)), mean(abs(diff_matrix(:))./slt_real(:), "omitmissing"));
fprintf('Absolute differences = |nfaslt_slt - |nfasltcx_slt|^2|\n')
fprintf('\tMin: %.2e, Max: %.2e, Mean: %.2e, Std: %.2e\n', ...
    min(abs(diff_matrix(:))), max(abs(diff_matrix(:))), mean(abs(diff_matrix(:))), std(abs(diff_matrix(:))));

% Check specific frequencies
fprintf('\nTarget frequency analysis:\n');
for f_target = vfTarget
    [~, idx] = min(abs(frequency_vector - f_target));
    actual_f = frequency_vector(idx);
    real_max = max(slt_real(idx, :));
    complex_max = max(slt_power_from_complex(idx, :));
    diff_max = max(abs(diff_matrix(idx, :)));
    fprintf('%.0f Hz (actual %.0f Hz): Real max=%.4f, Complex^2 max=%.4f, Max diff=%.2e\n', ...
        f_target, actual_f, real_max, complex_max, diff_max);
end

% Check for NaN or Inf values
real_nan = sum(~isfinite(slt_real(:)));
complex_nan = sum(~isfinite(slt_complex(:)));
fprintf('\nNon-finite values - Real: %d, Complex: %d\n', real_nan, complex_nan);

%% Machine precision tolerance

fprintf('\nValidation result:\n');
fprintf('Maximum absolute difference: %.2e\n', max_abs_diff);
fprintf('Maximum relative difference: %.2e\n', max_rel_diff);
fprintf('Tolerance threshold:         %.2e\n', tolerance);

if max_abs_diff < tolerance
    fprintf('✓ TEST PASSED: Functions agree within machine precision\n');
else
    fprintf('✗ TEST FAILED: Difference exceeds tolerance\n');
end

%% Visual comparison

time = linspace(0, N, numel(xSignal)) / 1000 - nPreSpaceS * fs / 1000;

figure;
tiledlayout(3, 2)

% Original time domain signal
nexttile(1)
plot(time, xSignal);
title('Test Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Target frequency power over time
nexttile(2)
colors = {'r', 'g', 'b'};
hold on;
for i = 1:length(vfTarget)
    [~, idx] = min(abs(frequency_vector - vfTarget(i)));
    plot(time, slt_real(idx, :), [colors{i} '-'], 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%.0f Hz (real)', vfTarget(i)));
    plot(time, slt_power_from_complex(idx, :), [colors{i} 'o--'], 'LineWidth', 1, ...
         'DisplayName', sprintf('%.0f Hz (|complex|^2)', vfTarget(i)));
end
legend('Location', 'best');
title('Target Frequencies Comparison');
xlabel('Time (s)');
ylabel('Power');
grid on;

% NFASLT Power Scalogram
nexttile(3)
imagesc(time, frequency_vector, slt_real);
set(gca, 'ydir', 'normal');
colormap(gca, jet);
title('nfaslt (real power)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% NFASLTCX Power Scalogram
nexttile(4)
imagesc(time, frequency_vector, slt_power_from_complex);
set(gca, 'ydir', 'normal');
colormap(gca, jet);
title('|nfasltcx|^2 (power from complex)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Power Scalogram Difference
nexttile(5)
imagesc(time, frequency_vector, diff_matrix);
set(gca, 'ydir', 'normal');
colormap(gca, jet);
title('Difference (real - |complex|^2)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% NFASLTCX Phase Scalogram
nexttile(6)
imagesc(time, frequency_vector, angle(slt_complex));
set(gca, 'ydir', 'normal');
colormap(gca, jet);
title('Phase from nfasltcx');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

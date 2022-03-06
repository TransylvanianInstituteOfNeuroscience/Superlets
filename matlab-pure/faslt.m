%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FRACTIONAL ADAPTIVE SUPERRESOLUTION WAVELET (SUPERLET) TRANSFORM 
% 
%   AUTHOR:         Harald Bârzan
%   DATE:           April 2019
%   DESCRIPTION:
%
%   Computes the adaptive superresolution wavelet (superlet) transform on 
%   input data to produce a time-frequency representation. For each 
%   frequency of interest, the closest integer order from the order 
%   interval will be chosen to produce each superlet. A superlet is a set 
%   of wavelets with the same center frequency but different number of 
%   cycles.
%
%   REFERENCE:
%   
%   Time-frequency super-resolution with superlets
%   Moca, V.V., Nagy-Dãbâcan, A., Bârzan, H., Mure?an, R.C.
%   https://www.nature.com/articles/s41467-020-20539-9
%   
%   NOTES:
%
%   If the input data consists of multiple buffers, a wavelet spectrum will
%   be computed for each of the buffers and averaged to produce the final 
%   result.
%   If the order parameter (ord) is empty, this function will return the
%   standard CWT (one wavelet per frequency of interest).
%
%   INPUT:
%   > input         - [buffers x samples] matrix
%   > Fs            - sampling frequency in Hz
%   > F             - frequency-of-interest buffer
%   > c1            - number of initial wavelet cycles
%   > o             - [1 x 2] interval of superresolution orders (optional)
%   > mult          - specifies the use of multiplicative superresolution
%                     (0 - additive, != 0 - multiplicative)
%
%   OUTPUT:
%   > wtresult      - [frequencies x samples] superlet spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wtresult] = faslt(input, Fs, F, c1, o, mult)

% check frequency of interest parameter
if (isempty(F))
    error('frequencies not defined'); 
end

% check order parameter and initialize the order used at each frequency. 
% if empty, go with an order of 1 for each frequency (single wavelet per
% set)
if (~isempty(o))
    % fractional order is whatever the user put in
    order_frac  = linspace(o(1), o(2), numel(F));
    % integer order is the real number of wavelets in the set
    order_int   = ceil(order_frac);
else
    order_frac  = ones(numel(F), 1);
    order_int   = order_frac;
end

% validate input buffer
if (isempty(input))
    error('input is empty'); 
end

% if input is a column vector, turn it into a row vector instead
if (size(input, 2) == 1 && size(input, 1) > 1)
    input = input'; 
end

% get the input size
[Nbuffers, Npoints] = size(input);

% the padding will be size of the lateral zero-pads, which serve to avoid
% border effects during convolution
padding = 0;

% the wavelet sets
wavelets = cell(numel(F), max(order_int));
        
% initialize wavelet sets for either additive or multiplicative
% superresolution

for i_freq = 1 : numel(F)
    for i_ord = 1 : order_int(i_freq)
        
        % compute the number of cycles (additive or multiplicative)
        if (mult ~= 0)
            n_cyc = i_ord * c1;
        else
            n_cyc = i_ord + c1;
        end
        
        % add the wavelet to the set
        wavelets{i_freq, i_ord} = cxmorlet(F(i_freq), n_cyc, Fs);
        
        % the margin will be the half-size of the largest wavelet
        padding = max(padding, fix(numel(wavelets{i_freq, i_ord}) / 2));
    end
end

% the zero-padded buffer
buffer = zeros(Npoints + 2 * padding, 1);

% the output scalogram
wtresult = zeros(numel(F), Npoints);

% convenience indexers for the zero-padded buffer
bufbegin    = padding + 1;
bufend      = padding + Npoints;

% loop over the input buffers
for i_buf = 1 : Nbuffers
    for i_freq = 1 : numel(F)
        % pooling buffer, starts with 1 because we're doing geometric mean
        temp = ones(1, Npoints);
        
        % fill the central part of the buffer with input data
        buffer(bufbegin : bufend) = input(i_buf, :);
        
        % get the number of integer wavelets
        n_wavelets = floor(order_frac(i_freq));
        
        % compute the convolution of the buffer with each wavelet in the
        % current set (integer wavelets)
        for i_ord = 1 : n_wavelets
            % restricted convolution (input size == output size)
            tempcx = conv(buffer, wavelets{i_freq, i_ord}, 'same');
            
            % accumulate the magnitude (times 2 to get the full spectral
            % energy), pool with exponent = 1
            temp = temp .* (2 .* abs(tempcx(bufbegin : bufend)) .^ 2)';
        end
        
        % handle fractional exponent
        if (is_fractional(order_frac(i_freq)) && ...
            ~isempty(wavelets{i_freq, order_int(i_freq)}))
            % set the order index
            i_ord = order_int(i_freq);
            
            % the exponent is the fractional remainder
            exponent = order_frac(i_freq) - fix(order_frac(i_freq));
            
             % restricted convolution (input size == output size)
            tempcx = conv(buffer, wavelets{i_freq, i_ord}, 'same');
            
            % accumulate the magnitude (times 2 to get the full spectral
            % energy), pool with exponent = 1
            temp = temp .* ((2 .* abs(tempcx(bufbegin : bufend)) .^ 2)') .^ exponent;
        end
            
        % compute the order of the geometric mean
        root = 1 / order_frac(i_freq);
        temp = temp .^ root;
        
        % accumulate the current FOI to the result spectrum
        wtresult(i_freq, :) = wtresult(i_freq, :) + temp;
    end
end

% scale the output by the number of input buffers
wtresult = wtresult ./ Nbuffers;

return


% computes the complex Morlet wavelet for the desired center frequency Fc
% with Nc cycles, with a sampling frequency Fs.
function w = cxmorlet(Fc, Nc, Fs)
    %we want to have the last peak at 2.5 SD
    sd  = (Nc / 2) * (1 / Fc) / 2.5;
    wl  = 2 * floor(fix(6 * sd * Fs)/2) + 1;
    w   = zeros(wl, 1);
    gi  = 0;
    off = fix(wl / 2);
    
    for i = 1 : wl
        t       = (i - 1 - off) / Fs;
        w(i)    = bw_cf(t, sd, Fc);
        gi      = gi + gauss(t, sd);
    end
    
    w = w ./ gi;
return

% compute the complex wavelet coefficients for the desired time point t,
% bandwidth bw and center frequency cf
function res = bw_cf(t, bw, cf)
    cnorm   = 1 / (bw * sqrt(2 * pi));
    exp1    = cnorm * exp(-(t^2) / (2 * bw^2));
    res     = exp(2i * pi * cf * t) * exp1;
return;

% compute the gaussian coefficient for the desired time point t and
% standard deviation sd
function res = gauss(t, sd)
    cnorm   = 1 / (sd * sqrt(2 * pi));
    res     = cnorm * exp(-(t^2) / (2 * sd^2));
return;

% tell me if a number is an integer or a fractional
function res = is_fractional(x)
    res = fix(x) ~= x;
return;
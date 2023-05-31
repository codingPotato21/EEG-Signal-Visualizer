%% Initial setup

close all
clear all
clc

load('EEG.mat');

%% Get the signals

% Get the left and right EEG signals from file
time = b1(:,1);
left_eeg = b1(:,2);
right_eeg = b1(:,3);

figure
plot(time, left_eeg)
title("Original signal")

%% Find the FFT

% Find the shifted FFT
n = length(left_eeg);
signal_fft = fftshift(abs(fft(left_eeg, n)));
frequencies = (0:n-1);

% Find the center peak of the FFT
[peaks, positions] = findpeaks(signal_fft, 'SortStr', 'descend');
max_positions = positions(1:2);
center = (max_positions(2) + max_positions(1)) / 2;

% Find the peaks where the noise happen
peakThresh = peaks(1) * 0.25;
[peaks, positions] = findpeaks(signal_fft, 'MinPeakHeight', peakThresh);
first_peak = min(positions);
last_peak = max(positions);

figure
plot(frequencies, signal_fft)
hold('on')
plot([first_peak last_peak max_positions(1) max_positions(2)], [0 0 0 0], 'rv', 'MarkerFaceColor', 'b')
hold('off')
title("Signal Peaks")

%% Center the FFT around 0

% Center the shifted FFT at 0 on the x-axis
range = (length(signal_fft) - center) * 2;
shifted_frequencies = range*linspace(-1, 1, n);

figure
plot(shifted_frequencies, signal_fft);
title("Shifted FFT")

%% Filter the noise out using bandpass filter

% Find sampling freq
fs = n / max(time);

% Filter parameters
shift_region = fs * 0.075;
fstop1 = (first_peak / 100) - (shift_region * 2);
fstop2 = (first_peak / 100) + (shift_region * 2);
fpass1 = fstop1 + (shift_region);
fpass2 = fstop2 - (shift_region);

% Find filter order
Rp = 1; % passband ripple
Rs = 40; % stopband attenuation
Wp = [fpass1, fpass2] / (fs/2); % Normalized passband freq
Ws = [fstop1, fstop2] / (fs/2); % Normalized stopband freq
[order, Wn] = buttord(Wp, Ws, Rp, Rs);

% Design the filter
[b, a] = butter(order, Wn, 'stop');

% Filter the original signal
filtered_signal = filter(b, a, left_eeg);

% Find fft of filtered signal
filtered_fft = fftshift(abs(fft(filtered_signal, n)));

figure
plot(shifted_frequencies, filtered_fft);
title("Filtered FFT")

%% Find the gamma signal using bandpass filter

% Get the filtered signal

% Filter parameters
fpass1 = 32;
fstop1 = 31.5;

% Find filter order
Rp = 1;
Rs = 60;
Wp = fpass1 / (fs/2); % Normalized passband freq
Ws = fstop1 / (fs/2); % Normalized stopband freq
[order, Wn] = cheb2ord(Wp, Ws, Rp, Rs, 's');

% Find filter coefficients
[b, a] = cheby2(order, Rs, Wn, 'high', 's');

filtered_signal = filter(b, a, filtered_fft);

figure
plot(time, filtered_signal)
title("Extracted gamma signal")

%% Find the alpha

% Filter parameters
% The frequencies for gamma signal is usually between 30 - 100 Hz
fpass1 = 8;
fpass2 = 12;
fstop1 = 7.5;
fstop2 = 12.5;

% Find filter order
Rp = 1;
Rs = 20;
Wp = [fpass1, fpass2] / (fs/2); % Normalized passband freq
Ws = [fstop1, fstop2] / (fs/2); % Normalized cutoff freq
[order, Wn] = cheb2ord(Wp, Ws, Rp, Rs, 's');

% Find filter coefficients
[b, a] = cheby2(order, Rs, Wn, 'bandpass', 's');

filtered_signal = filter(b, a, filtered_fft);

figure
plot(time, filtered_signal)
title("Extracted alpha signal")

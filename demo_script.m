%% Initial setup

close all
clear all
clc

load('EEG.mat');

%% Get the signals

% Get the left and right EEG signals from file
file_time = b1(:,1);
left_eeg = b1(:,2);
right_eeg = b1(:,3);

% Trim the signal and remove 0 values at the start and end
first_nonzero_index = find(left_eeg ~= 0, 1, 'first');
last_nonzero_index = find(left_eeg ~= 0, 1, 'last');
left_eeg = left_eeg(first_nonzero_index:last_nonzero_index);

% Find the sampling freq and signal time
n = length(left_eeg);
fs = n / max(file_time);
time = (0:n-1) / fs;

fg = figure();
subplot(4, 1, 1)
plot(time, left_eeg)
title("Original signal")
ylabel("Voltage (uV)")
xlabel("Time (seconds)")

%% Find the FFT

n = length(left_eeg);
left_fft = abs(fft(left_eeg, n));
frequencies = (0:n-1) * (fs/n);

subplot(4, 1, 2)
plot(frequencies, left_fft)
title("Original FFT")
ylabel("|x|")
xlabel("Frequency (Hz)")

%% Filter the original signal using frequencies from the FFT

% Filter signal with filter designer

Fpass1 = 48.5;        % First Passband Frequency
Fstop1 = 49;          % First Stopband Frequency
Fstop2 = 51;          % Second Stopband Frequency
Fpass2 = 51.5;        % Second Passband Frequency
Apass1 = 1;           % First Passband Ripple (dB)
Astop  = 60;          % Stopband Attenuation (dB)
Apass2 = 1;           % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

filtered_signal = filter(Hd, left_eeg);

filtered_fft = abs(fft(filtered_signal));

subplot(4, 1, 3)
plot(time, filtered_signal);
title("Filtered Signal")
ylabel("Voltage (uV)")
xlabel("Time (seconds)")

subplot(4, 1, 4)
plot(frequencies, filtered_fft);
title("Filtered FFT")
ylabel("|x|")
xlabel("Frequency (Hz)")

%% Gamma signal

Fstop = 31.5;        % Stopband Frequency
Fpass = 32;          % Passband Frequency
Astop = 80;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

gamma_signal = filter(Hd, filtered_signal);

fg2 = figure();
subplot(5, 1, 1)
plot(time, gamma_signal)
title("Gamma Signal")
xlabel("Time (seconds)")
ylabel("Voltage (uV)")
xlim([0, 10])
ylim('auto')

%% Beta signal

Fstop1 = 12.5;        % First Stopband Frequency
Fpass1 = 13;          % First Passband Frequency
Fpass2 = 30;          % Second Passband Frequency
Fstop2 = 30.5;        % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

beta_signal = filter(Hd, filtered_signal);

subplot(5, 1, 2)
plot(time, beta_signal)
title("Beta Signal")
xlabel("Time (seconds)")
ylabel("Voltage (uV)")
xlim([0, 10])
ylim('auto')

%% Alpha signal

Fstop1 = 7.5;         % First Stopband Frequency
Fpass1 = 8;           % First Passband Frequency
Fpass2 = 12;          % Second Passband Frequency
Fstop2 = 12.5;        % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

alpha_signal = filter(Hd, filtered_signal);

subplot(5, 1, 3)
plot(time, alpha_signal)
title("Alpha Signal")
xlabel("Time (seconds)")
ylabel("Voltage (uV)")
xlim([0, 10])
ylim('auto')

%% Theta signal

Fstop1 = 3.5;         % First Stopband Frequency
Fpass1 = 4;           % First Passband Frequency
Fpass2 = 7;           % Second Passband Frequency
Fstop2 = 7.1;         % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

theta_signal = filter(Hd, filtered_signal);

subplot(5, 1, 4)
plot(time, theta_signal)
title("Theta Signal")
xlabel("Time (seconds)")
ylabel("Voltage (uV)")
xlim([0, 10])
ylim('auto')

%% Delta signal

Fstop1 = 0.01;        % First Stopband Frequency
Fpass1 = 0.2;         % First Passband Frequency
Fpass2 = 4;           % Second Passband Frequency
Fstop2 = 4.1;         % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

delta_signal = filter(Hd, filtered_signal);

subplot(5, 1, 5)
plot(time, delta_signal)
title("Delta Signal")
xlabel("Time (seconds)")
ylabel("Voltage (uV)")
xlim([0, 10])
ylim('auto')

%% Plotting filter

% Compute the frequency response
freq = linspace(0, fs/2, n);
[H, f] = freqz(Hd, freq, fs);

% Convert freq response to phase response in radians
phase = unwrap(angle(H));

% Find the maximum value of the phase response
max_phase = max(phase);

% Align the phase response by subtracting the maximum phase value
aligned_phase = phase - max_phase;

% Plot magnitude response in dB
figure
yyaxis left;
plot(f, 20*log10(H), 'b');
xlabel('Frequency (Hz)');
ylabel('Magnitude Response (dB)');
title('Delta Filter Magnitude/Phase Response');

% Set the y-axis limits for each plot
ylim(gca, [-100, 0]); % Set the y-axis limits for the magnitude plot

yyaxis right
% Plot phase response in radians
plot(f, aligned_phase, 'r')
ylabel("Phase Response (radians)")

% Set the y-axis limits for each plot
ylim(gca, [min(aligned_phase), 0]); % Set the y-axis limits for the magnitude plot

xlim(gca, [0, fs/2])
grid("on")
legend("Magnitude", "Phase")

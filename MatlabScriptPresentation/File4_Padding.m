% script by STO ten Oever, 09-02-2018. Developed for M-BIC: Disclosing fine-grained temporal 
% processing: Common and advanced analysis of EEG time-series. 
% 
% script 4. 
% This script provides examples of padding effects
%

%% example data with two frequencyies and padded or unpadded data
freqD = 5;
freqD2 = 8; % !!! try different frequencies here to see effect.
fsample = 100;
tp = [1/fsample:1/fsample:1.1];
ndatsample = length(tp);
ampad = ndatsample*2; % !!! try various padding to see effect.
ndatsampleP = ndatsample + ampad*2;

dat = sin(freqD*2*pi*tp)+sin(freqD2*2*pi*tp);
dat = dat - mean(dat);
dat_pad = [zeros(1,ampad) dat zeros(1,ampad)];
freq = linspace(0, (ndatsample/2)./(ndatsample./fsample), ndatsample/2+1);
freq_p = linspace(0, (ndatsampleP/2)./(ndatsampleP./fsample), ndatsampleP/2+1);

tap = hanning(ndatsample)';
tap = tap./norm(tap, 'fro');

tap_p = hanning(ndatsampleP)';
tap_p = tap_p./norm(tap_p, 'fro');

%% Show the effect on the FFT 
close all
figure(1)
subplot(121)
fft_output = fft(bsxfun(@times,dat,tap),[], 2);
fft_output = fft_output ./norm(fft_output , 'fro');
plot(freq, abs(fft_output(1:ndatsample/2+1)).^2, '.-', 'markersize', 10); xlabel('Freq (Hz)'); ylabel('amplitude (au)');
title('normal'); set(gca, 'xlim', [0 10], 'ylim', [0 0.3]);
hold on
plot([freqD freqD], [0 0.5]);
plot([freqD2 freqD2], [0 0.5]);
subplot(122)
fft_output_pad = fft(bsxfun(@times,dat_pad,tap_p),[], 2);
fft_output_pad = fft_output_pad ./norm(fft_output_pad , 'fro');
plot(freq_p, abs(fft_output_pad(1:ndatsampleP/2+1)).^2, '.-', 'markersize', 10); xlabel('Freq (Hz)'); ylabel('amplitude (au)');
title('padded'); set(gca, 'xlim', [0 10], 'ylim', [0 0.3]);
hold on
plot([freqD freqD], [0 0.5]);
plot([freqD2 freqD2], [0 0.5]);
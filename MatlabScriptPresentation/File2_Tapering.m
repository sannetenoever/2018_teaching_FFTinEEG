% script by STO ten Oever, 09-02-2018. Developed for M-BIC: Disclosing fine-grained temporal 
% processing: Common and advanced analysis of EEG time-series. 
% 
% script 2. 
% This script provides examples of tapering
%
% necessary before usage:
%       - add fieldtrip to your path (change line 10).

addpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20180712');

%% Create data and put also in Fieldtrip structure
fsample = 100;
tp = [1/fsample:1/fsample:0.9];
freq = 5;
ndatsample = length(tp);

dat = sin(freq*2*pi*tp)+rand(1,ndatsample)*0;
dat = dat - mean(dat);

% make fake fieldtrip structure
fftr.trial{1} = dat;
fftr.time{1} = tp;
fftr.label{1} = 'fake';

figure(1)
subplot(231)
plot(tp, dat); xlabel('time (sec)'); ylabel('amplitude (au)'); set(gca, 'ylim', [-1.1 1.1]);
title('main data');

%% perform FFT with Fieldtrip
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'pow';
fr = ft_freqanalysis(cfg, fftr);

figure(1)
subplot(232)
plot(fr.freq, sqrt(fr.powspctrm), '.-', 'markersize', 20); xlabel('freq (Hz)'); ylabel('amplitude (au)');
title('fieldtrip output');
% amplitude is shown instead of power [ sqrt(power) = amplitude ] as effect is clearer

%% perform fourier analysis with matlab fft
fft_output = fft(dat,[], 2);
fft_output = fft_output ./norm(fft_output, 'fro');
fft_output = abs(fft_output);

figure(1)
subplot(233)
plot(fr.freq, fft_output(1:ceil(ndatsample/2)+1), '.-', 'markersize', 20); xlabel('freq (Hz)'); ylabel('amplitude (au)');
title('FFT output');

%% make a hanning taper
tap = hanning(ndatsample)';

figure(1)
subplot(234)
plot(tp, tap); xlabel('freq (Hz)'); ylabel('power (au)');
title('taper');

figure(1)
subplot(235)
plot(tp, bsxfun(@times,dat,tap)); xlabel('time (sec)'); ylabel('amplitude (au)'); set(gca, 'ylim', [-1.1 1.1]);

%% perform fourier analysis with matlab fft using the hanning tapered data
fft_output_tap = fft(bsxfun(@times,dat,tap),[], 2);
fft_output_tap = fft_output_tap ./norm(fft_output_tap , 'fro');
fft_output_tap = abs(fft_output_tap);

figure(1)
subplot(236)
plot(fr.freq, fft_output_tap(1:ceil(ndatsample/2)+1), '.-', 'markersize', 20); xlabel('freq (Hz)'); ylabel('amplitude (au)');
title('FFT output'); 


%% Explain periodic assumption and why tapering is necessary
freqs = [5 5.556];

for fr = 1:2
    fsample = 100;
    tp = [1/fsample:1/fsample:0.9];
    ndatsample = length(tp);
    dat = sin(freqs(fr)*2*pi*tp)+rand(1,ndatsample)*0;
    dat = dat - mean(dat);
    
    fft_output = fft(dat,[], 2);
    fft_output = fft_output ./norm(fft_output, 'fro');
    fft_output = abs(fft_output); % amplitude
    
    % go to the frequency we have in the FFT spectrum
    freqboilim = round([0 ndatsample/2]) + 1; % limits in fourierpoints, determined by the window length.
    freqboi    = freqboilim(1):1:freqboilim(2);
    freqoi     = (freqboi-1) ./ (ndatsample./fsample); % frequency is determined by sampling rate
    
    figure(6)
    subplot(2,2, fr);
    plot(freqoi, fft_output(1:ceil(ndatsample/2)+1), '.-', 'markersize', 20); xlabel('freq (Hz)'); ylabel('amplitude (au)');
    title('FFT output');
    subplot(2,2, fr+2)
    plot([tp tp+tp(end)], [dat dat]); xlabel('time (sec)'); ylabel('amplitude (au)');
    set(gca, 'ylim', [-1.1 1.1]);
end;

%% Provide example with different tapers in amplitude + decibel amplitude
fsample = 100;
tp = [1/fsample:1/fsample:6.9];
freq = 5;
ndatsample = length(tp);

dat = sin(freq*2*pi*tp)+rand(1,ndatsample)*0;
dat = dat - mean(dat);
FreqUse = [0:ndatsample/2]./tp(end); 

% different types of tapers
tapd{1} = ones(1,ndatsample);
tapd{2} = hanning(ndatsample)';
tapd{3} = [zeros(1,ceil(ndatsample/4)) ones(1,floor(ndatsample/4)) ones(1,floor(ndatsample/4)) zeros(1,ceil(ndatsample/4))];
tapd{4} = cos(linspace(-1, 1,ndatsample));
tapd{5} = tukeywin(ndatsample,0.5)';
tapd{6} = hamming(ndatsample)';

names = {'none';'hanning';'square';'cos';'tukeywin0.5';'hamming'};
clear saveoutd saveout

% calculate the output:
for it = 1:length(tapd)
    tapd{it} = tapd{it}./norm(tapd{it}, 'fro');
    fft_output_tap = fft(bsxfun(@times,dat,tapd{it}),[], 2);
    saveoutd(it,:) = fft_output_tap;
    fft_output_tap = fft_output_tap ./norm(fft_output_tap , 'fro');
    fft_output_tap = abs(fft_output_tap);
    saveout(it,:) = fft_output_tap;
end;

% plot amplitude
for it = 1:length(tapd)
    figure(2)
    subplot(6,2,1+(it-1)*2)
    plot(tp, tapd{it}); xlabel('time (sec)');set(gca, 'xlim', [0 6.9]);
    if it == length(tapd) 
        xlabel('time (sec)'); ylabel('amplitude (au)');title(names{it});
    end;
    subplot(6,2,it*2)   
    hold on
    plot(FreqUse, saveout(it,1:ceil(ndatsample/2)+1), '.-', 'markersize', 20); 
    if it == length(tapd)
        xlabel('freq (Hz)'); ylabel('amplitude (au)');
    end;
    title(names{it});
    set(gca, 'xlim', [4 8], 'ylim', [0 0.06]);
end;

% plot dB of amplitude
for it = 1:length(tapd)
    figure(3)
    subplot(6,2,1+(it-1)*2)
    plot(tp, tapd{it}); 
    if it == length(tapd) 
        xlabel('time (sec)'); ylabel('amplitude (au)');title(names{it});
    end;
    set(gca, 'xlim', [0 6.9]);
    subplot(6,2,it*2)   
    hold on
    plot(FreqUse, mag2db(saveout(it,1:ceil(ndatsample/2)+1)), '.-', 'markersize', 20); 
    if it == length(tapd)
        xlabel('freq (Hz)'); ylabel('dB amplitude (au)');
    end;
    title(names{it});
    set(gca, 'xlim', [4 8], 'ylim', [-100 0]);
end;

%% DPSS multitaper
fsample = 100;
tp = [1/fsample:1/fsample:3.9];
freq = 5;
ndatsample = length(tp);

dat = sin(freq*2*pi*tp)+rand(1,ndatsample)*0; % you can make the data noisy by changing 0 to a higher amplitude (such as 5)
dat = dat - mean(dat);
FreqUse = [0:ndatsample/2]./tp(end); 

% this is what you do with Fieldtrip to get the dpss multitaper:
% cfg = [];
% cfg.method = 'mtmfft';
% cfg.taper = 'dpss';
% cfg.tapsmofrq = 1;
% cfg.output = 'pow';
% fr = ft_freqanalysis(cfg, fftr);

% dpss taper: amount of tapers dependent on the amount of smoothing!
tap = dpss(ndatsample, double(ndatsample).*(cfg.tapsmofrq./fsample))';
% remove the last taper because the last slepian taper is always messy
tap = tap(1:(end-1), :);

% plot the taper, tapered data and FFT output
figure(4)
for it = 1:size(tap,1)
    subplot(size(tap,1),3,(it-1).*3+1);
    plot(tp, tap(it,:));
    if it == size(tap,1);
        xlabel('time (sec)'); ylabel('amplitude');
    end;
    subplot(size(tap,1),3,(it-1).*3+2);   
    plot(tp, bsxfun(@times,dat,tap(it,:)));    
    alltapdat(it,:) = bsxfun(@times,dat,tap(it,:));
    if it == size(tap,1);
        xlabel('time (sec)');
    end;  
    subplot(size(tap,1),3,(it-1).*3+3);    
    fft_output_tap = fft(bsxfun(@times,dat,tap(it,:)),[], 2);
    fft_output_tap = fft_output_tap ./norm(fft_output_tap , 'fro');
    fft_output_tap = abs(fft_output_tap).^2;
    allpow(it,:) = fft_output_tap(1:ceil(ndatsample/2)+1);
    plot(FreqUse, fft_output_tap(1:ceil(ndatsample/2)+1)); 
    set(gca, 'xlim', [0 15]);
    if it == size(tap,1);
        xlabel('Freq (Hertz)'); 
    end;  
end;

%% Plot tapers for different smoothing factor (=amount of tapers)
% first for pure sine, then for noisy data
close(figure(5))

datr = sin(freq*2*pi*tp);
allsm = [1 2 5]; % all smoothing factors
for it2 = 1:2;
    figure(5)
    subplot(1,2,it2)
    for it = 1:length(allsm)
        tap = dpss(ndatsample, double(ndatsample).*(allsm(it)./fsample))';
        tap = tap(1:(end-1), :);
        fft_output_tap = fft(bsxfun(@times,datr, tap),[], 2);
        fft_output_tap = fft_output_tap ./norm(fft_output_tap , 'fro');
        fft_output_tap = abs(fft_output_tap).^2;
        plot(FreqUse, sqrt(mean(fft_output_tap(:,1:ceil(ndatsample/2)+1))));
        xlabel('Freq (Hz)'); ylabel('amplitude (au)'); set(gca, 'xlim', [0 50]);
        hold on
    end;
    tap = hanning(ndatsample)';
    fft_output_tap = fft(bsxfun(@times,datr, tap),[], 2);
    fft_output_tap = fft_output_tap ./norm(fft_output_tap , 'fro');
    fft_output_tap = abs(fft_output_tap).^2;
    plot(FreqUse, sqrt(fft_output_tap(:,1:ceil(ndatsample/2)+1)));
    
    datr = datr+(rand(size(datr))-0.5).*5; % add noise
end;
leg = {'1 Hz smooth';'2 Hz smooth'; '5 Hz smooth'; 'hanning'};
legend(leg);



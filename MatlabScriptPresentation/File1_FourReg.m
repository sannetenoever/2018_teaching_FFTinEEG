% script by STO ten Oever, 09-02-2018. Developed for M-BIC: Disclosing fine-grained temporal 
% processing: Common and advanced analysis of EEG time-series. 
% 
% script 1. 
% This script provides an explantion of FFT by using sine and cosine fitting.
%

%% Create data
fsample = 100;
tp = [1/fsample:1/fsample:1.2];
freq = 5;
ndatsample = length(tp);

dat = sin(freq*2*pi*tp+0*pi); % add different pi values (phases) to see effect on sine and cosine fit
dat = dat - mean(dat);

figure(1)
plot(tp, dat); xlabel('time (sec)'); ylabel('amplitude (au)'); set(gca, 'ylim', [-2.1 2.1]);
title('main data');
 
%% calculate fft
fft_output = fft(dat,[], 2);
fft_outputP = abs(fft_output).^2; % pow
fft_output = abs(fft_output); % amplitude

figure(2)
subplot(131)
plot(fft_output, '.-', 'markersize', 20); ylabel('amplitude (au)');
title('FFT output');

set(gca,'xlim', [0 ndatsample], 'ylim', [-70 70]);

%% multiple data with sine for all integer datapoints 
clear RM bs f cov
ADP = linspace(1/ndatsample,ndatsample, ndatsample);
for it = 0:ndatsample % go through different 'frequencies'
    bh(it+1) = sum(dat.*sin(it*2*pi*ADP)); 
end;

figure(2)
subplot(132)
plot(bh, '.-', 'markersize', 20); ylabel('amplitude (au)');
title('SinusFit');

set(gca,'xlim', [0 ndatsample], 'ylim', [-70 70]);

% show data, sine example (it = 5) and their multiplication
figure(3)
it = 5;
subplot(211)
plot(dat); hold on 
plot(sin(it*2*pi*ADP));
legend({'data'; ['sine: it ' num2str(it)]});
subplot(212)
plot(dat.*sin(it*2*pi*ADP));

%% multiple data with cosine for all integer datapoints 
clear RM bs f
for it = 0:ndatsample % go through different 'frequencies'
    bh(it+1) = sum(dat.*cos(it*2*pi*ADP)); 
end;

figure(2)
subplot(133)
plot(bh, '.-', 'markersize', 20); ylabel('amplitude (au)');
title('CosFit');

set(gca,'xlim', [0 ndatsample], 'ylim', [-70 70]);

%% multiply data with both sine and cosine for all integer datapoints 
clear RM bs f
for it = 0:ndatsample
    ah(it+1) = sum(dat.*cos(it*2*pi*ADP));
    bh(it+1) = sum(dat.*sin(it*2*pi*ADP));    
    f(it+1) = ah(it+1)+i*bh(it+1);
end;

figure(3)
subplot(121)
plot(fft_output, '.-', 'markersize', 20); ylabel('amplitude (au)');
title('FFT output');
set(gca,'xlim', [0 ndatsample], 'ylim', [0 70]);

% plot the amplitude (abs(complexFourier))
subplot(122)
plot(abs(f), '.-', 'markersize', 20); ylabel('amplitude (au)');
title('sinus+cosine output');
set(gca,'xlim', [0 ndatsample], 'ylim', [0 70]);

%% Show the fitted sine for various values of it:
% the interesting effect occurs for very high it (frequencies)
itv = [1 2 3 4 6 15 40 ndatsample-4];
cntpl = 1;
for it = 1:length(itv)
    figure(4)
    subplot(2,length(itv)/2,cntpl);  
    plot(sin(itv(it).*2*pi*ADP)); title(['it=' num2str(itv(it)-1) '; sine'])
    set(gca, 'xlim', [0 ndatsample]);
    cntpl = cntpl+1;
end;

%% "aliasing" in the FFT, sin and cos
close(figure(5))
figure(5)
it = 115;
it2 = 5;
ADPo = linspace(1/ndatsample,ndatsample, ndatsample);
ADPw = linspace(1/ndatsample,ndatsample.*100, ndatsample.*100);
LSO = linspace(1/ndatsample,ndatsample, ndatsample.*100);
So =  sin(it.*2*pi*ADPo);
So2 = sin(it2.*2*pi*ADPo);
Sw = sin(it.*2*pi*ADPw);
subplot(2,3,1);
plot(LSO, Sw);    
hold on
plot(ADPo, So, '.-','markersize', 15);
set(gca, 'xlim', [0 4]);
subplot(2,3,2);
plot(LSO, Sw);    
hold on
plot(ADPo, So, '.-','markersize', 15); title('sine');
set(gca, 'xlim', [1/ndatsample, ndatsample/4]);
subplot(2,3,3);
plot(ADPo, So2, '.-','markersize', 15);
set(gca, 'xlim', [1/ndatsample, ndatsample/4]);
Co =  cos(it.*2*pi*ADPo);
Co2 = cos(it2.*2*pi*ADPo);
Cw = cos(it.*2*pi*ADPw);
subplot(2,3,4);
plot(LSO, Cw);    
hold on
plot(ADPo, Co, '.-','markersize', 15);
set(gca, 'xlim', [0 4]);
subplot(2,3,5);
plot(LSO, Cw);    
hold on
plot(ADPo, Co, '.-','markersize', 15); title('cosine');
set(gca, 'xlim', [1/ndatsample, ndatsample/4]);
subplot(2,3,6);
plot(ADPo, Co2, '.-','markersize', 15);
set(gca, 'xlim', [1/ndatsample, ndatsample/4]);

%% calculation for frequencies that we estimated
amountCy = 0:ndatsample; % amount of cycles we try to put in the data window
FreqUse = amountCy./tp(end); % length of data window > go to Hertz
FreqUse = FreqUse(1:ndatsample/2+1); % only half is needed
FreqRes = FreqUse(2)-FreqUse(1); % frequency resolution
FreqRes = 1/tp(end); % can be directly calculated

% script by STO ten Oever, 09-02-2018. Developed for M-BIC: Disclosing fine-grained temporal 
% processing: Common and advanced analysis of EEG time-series. 
% 
% script 5. 
% This script provides explanation of phase vs power and weird fieldtrip
% angle behavior
%
% necessary before usage:
%       - add fieldtrip to your path (change line 10).

addpath('F:\Other\Teaching\MBicEEGcourse\Practical\MatlabToolboxes\fieldtrip-20161231');

%% Example data
fsample = 100;
tp = [1/fsample:1/fsample:1.1];
freq = 5;
ndatsample = length(tp);
FreqUse = linspace(0, (ndatsample/2)./(ndatsample./fsample), ndatsample/2+1);

dat = sin(freq*2*pi*tp)+rand(1,ndatsample).*0.1;
dat = dat - mean(dat);

figure(1)
plot(tp, dat); xlabel('time (sec)'); ylabel('amplitude (au)'); set(gca, 'ylim', [-2.1 2.1]);
title('main data');

%% calculate FFT, show power and phase
tap = hanning(ndatsample)';
fft_output = fft(bsxfun(@times,dat,tap),[], 2);
fft_output = fft_output ./norm(fft_output, 'fro');
fft_outputP = abs(fft_output).^2; % pow
fft_outputA = angle(fft_output);

figure(2)
subplot(121)
plot(fft_outputP, '.-', 'markersize', 20); ylabel('amplitude (au)');
title('pow');
subplot(122)
plot(fft_outputA, '.-', 'markersize', 20); ylabel('phase (rad)');
title('phase');

%% Show effect on phase estimate with different noise levels
noiselevel = [0 5 22];
amtr = 100;
intf = nearest(FreqUse, freq);

close(figure(3));
figure(3); cntpl = 1;
for nl = 1:length(noiselevel)
    Pow=zeros(1,ndatsample);PE5=[];
    for it = 1:amtr
        dat = sin(freq*2*pi*(tp)+0.5*pi)+rand(1,ndatsample).*noiselevel(nl);
        dat = dat - mean(dat);
        
        fft_output = fft(bsxfun(@times,dat,tap),[], 2);
        fft_output = fft_output ./norm(fft_output, 'fro');
        fft_outputP = abs(fft_output).^2; % pow
        fft_outputA = angle(fft_output);
        PE5(it) = fft_outputA(intf);
        Pow = Pow+fft_outputP;
    end;
    subplot(3,2,cntpl);
    plot(Pow./amtr, '.-');
    set(gca, 'xlim', [1 ndatsample./2]); xlabel('Freq'); ylabel('pow');
    subplot(3,2,cntpl+1);
    hist(PE5);
    set(gca, 'xlim', [-pi pi]); xlabel('phase (rad)'); ylabel('n');
    cntpl = cntpl + 2;
end;

%%
fsample = 1000;
freq = 5;
tp = [0:1/fsample:1-1/fsample];

ndatsample = length(tp);
FreqUse = linspace(0, (ndatsample/2)./(ndatsample./fsample), ndatsample/2+1);
intf = nearest(FreqUse, freq);
tap = hanning(ndatsample)';

dat = sin(freq*2*pi*tp+1.5*pi);
dat = dat - mean(dat);

fft_output = fft(bsxfun(@times,dat,tap),[], 2);
fft_outputP = abs(fft_output).^2; % pow
fft_outputA = angle(fft_output);
tp = tp-0; % !!!! change this to see weird fieltrip behavior, optimal is to change within one period, for 5 Hz between 0 and 0.2
fftr.time{1}= tp;
fftr.trial{1}= dat;
fftr.label{1} = '1';

cfg =[];
cfg.taper = 'hanning';
cfg.method = 'mtmfft';
cfg.output = 'fourier';
x = ft_freqanalysis(cfg, fftr);
FTph = angle(x.fourierspctrm(intf));
FFTph = fft_outputA(intf);

figure(3)
plot(fftr.time{1}, dat);
title(['FT = ' num2str(round(FTph./pi,2))  ' pi; FFT: ' num2str(round(FFTph./pi,2)) ' pi']);
set(gca, 'xlim', fftr.time{1}([1 end]));
xlabel('time (sec)');
hold on
plot([mean(tp) mean(tp)], [-1 1]);
plot([0 0], [-1 1]);
hold off

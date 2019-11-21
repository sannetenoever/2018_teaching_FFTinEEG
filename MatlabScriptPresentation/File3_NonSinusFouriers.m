% script by STO ten Oever, 09-02-2018. Developed for M-BIC: Disclosing fine-grained temporal 
% processing: Common and advanced analysis of EEG time-series. 
% 
% script 3. 
% This script provides examples of different kinds of data and their
% Fourier transform
%

clear

%% sample example data
%close all
fsample = 100;
tp = [1/fsample:1/fsample:2];
ndatsample = length(tp);
FreqUse = [0:ndatsample/2]./tp(end); 

%% non-sinus ones
clear dat
% linear
dat{1} = linspace(-1,1,ndatsample);
% impulse1
dat{2} = zeros(1,ndatsample);
dat{2}(round(length(tp)/2)) = 1;
dat{2}(round(length(tp)/2)+1) = -1;
% impulse2
dat{3} = zeros(1,ndatsample);
dat{3}(round(length(tp)/2)) = 1;
dat{3}(round(length(tp)/2)+1) = 1;
% sawtooth
dat{4} = sawtooth(tp*2*pi*5);
% sawtooth + linear trend
dat{5} = sawtooth(tp*2*pi*5)+dat{1}.*2;

%% non periodic sinus
clear dat
% 180 phase reset:
dat{1} = [sin(2*pi*5*tp(1:length(tp)/2)) sin(2*pi*5*tp(1:length(tp)/2)+pi)];
% 180 phase reset + power increase:
dat{2} = [sin(2*pi*5*tp(1:length(tp)/2)) sin(2*pi*5*tp(1:length(tp)/2)+pi).*2];
% increase in frequency:
t=linspace(0,10,length(tp));
dat{3} = [sin(2*pi/10.*t.^2)];
% amplitude modulation
y1=sin(2*pi*2*tp); % message signal
y2=sin(2*pi*20*tp); % carrier signal
dat{4} = (1+2.*y1).*(2.*y2);

%% untapered:
for it = 1:length(dat)
    dat{it} = dat{it}-mean(dat{it});
    fft_output = fft(dat{it},[], 2);
    fft_output = fft_output ./norm(fft_output, 'fro');
    fft_outputP = abs(fft_output).^2; % pow
    fft_output = abs(fft_output); % amplitude
        
    figure(1)
    subplot(length(dat),2,1+(it-1)*2)
    plot(tp, dat{it}); set(gca, 'xlim', [0 tp(end)]);
    if it == length(dat)
        xlabel('time (sec)'); ylabel('amplitude (au)');
    end;
    subplot(length(dat),2,it*2)
    hold on
    plot(FreqUse, fft_output(1:ceil(ndatsample/2)+1), '.-', 'markersize', 20);
    if it == length(dat)
        xlabel('freq (Hz)'); ylabel('amplitude (au)');
    end;
    hold off
    set(gca, 'xlim', [0 30]); 
end;

%% hanning tapered:
tap = hanning(ndatsample)';
tap = tap./norm(tap, 'fro');

for it = 1:length(dat)
    dat{it} = dat{it}-mean(dat{it});
    fft_output = fft(bsxfun(@times,dat{it},tap),[], 2);
    fft_output = fft_output ./norm(fft_output, 'fro');
    fft_outputP = abs(fft_output).^2; % pow
    fft_output = abs(fft_output); % amplitude
        
    figure(2)
    subplot(length(dat),2,1+(it-1)*2)
    plot(tp, dat{it}); set(gca, 'xlim', [0 tp(end)]);
    if it == length(dat)
        xlabel('time (sec)'); ylabel('amplitude (au)');
    end;
    subplot(length(dat),2,it*2)
    hold on
    plot(FreqUse, fft_output(1:ceil(ndatsample/2)+1), '.-', 'markersize', 20);
    if it == length(dat)
        xlabel('freq (Hz)'); ylabel('amplitude (au)');
    end;
    hold off
    set(gca, 'xlim', [0 30]);        
end;

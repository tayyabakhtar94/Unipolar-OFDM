%% inverse discrete fourier transform

% create the signal
srate = 1000; % hz
time = 0:1/srate:2; % time vector in seconds
pnts = length(time); % number of time points
signal = 2.5 * sin(2*pi*4*time) + 1.5 * sin (2*pi*6.5*time);

% prepare the fourier transform
fourTime = (0:pnts-1)/pnts;
fCoefs = zeros(size(signal));

for fi=1:pnts

% create complex sine wave
csw = exp(-1i*2*pi*(fi-1)*fourTime);

% compute dot product between sine wave and signal
fCoefs(fi) = sum(signal.*csw);
% these are called the fourier coefficients
end 

% extract amplitudes
ampls = abs(fCoefs) / pnts;
ampls(2:end) = 2*ampls(2:end);

% compute frequencies vector
hz = linspace(0,srate/2,floor(pnts/2)+1);

figure(1), clf
plot(hz,ampls(1:length(hz)),'s-')
% better:
stem(hz,ampls(1:length(hz)),'ks-','linew',3,'markersize',10,'markerfacecolor','w')

% make plot look a bit nicer
set(gca,'xlim',[0 10],'ylim',[-.01 3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u)')



%% inverse fourier transform

% initialize time-domain reconstruction 
reconsignal = zeros(size(signal));

for fi=1:pnts

% create coefficient-modulated complex sine wave
csw = fCoefs(fi) * exp(1i*2*pi*(fi-1)*fourTime);

% sum them together
reconsignal = reconsignal +csw;
end

% divide by N
reconsignal = reconsignal/pnts;

figure(2), clf
plot(time,signal)
hold on
plot(time,real(reconsignal),'ro')
legend({'original';'reconstructed'})
xlabel('Time (s)')
%% one frequency at a time 

clear

% set parametres
srate = 1000;
time = 0:1/srate:3;
pnts = length(time);

% create multispectral signal
signal = (1+sin(2*pi*12*time)) .* cos(sin(2*pi*25*time)+time);

% prepare the fourier transform
fourTime = (0:pnts-1)/pnts;
fCoefs = zeros(size(signal));

% here is the fourier transform
for fi=1:pnts
csw = exp(-1i*2*pi*(fi-1)*fourTime);
fCoefs(fi) = sum(signal.*csw)/pnts;
end

% frequencies in hz
hz = linspace(0,srate,pnts);

% setup plotting
figure(3), clf
subplot(211)
plot(time,signal), hold on
sigh = plot(time,signal,'k');
xlabel('Time (s)')
title('Time domain')

subplot(212)
powh = plot(1,'k','linew',2);
set(gca,'xlim',hz([1 end]),'xtick',0:100:900,'xticklabel',[0:100:500 400:-100]);
title('Frequency domain')
xlabel('Frequency (Hz)')

% initialize the reconstructed signal
reconsignal = zeros(size(signal));

% inverse fourier transform here
for fi=1:pnts

% create coeefficients-modulated complex sine wave
csw = fCoefs(fi) * exp(1i*2*pi*(fi-1)*fourTime);

% sum them together
reconsignal = reconsignal + csw;

% update plot for some frequencies
if  fi < dsearchn(hz',100) || fi>dsearchn(hz',srate-100)
set(sigh,'YData',real(reconsignal))
set(powh,'XData',hz(1:fi),'YData',2*abs(fCoefs(1:fi)))
pause(.05)
end
end
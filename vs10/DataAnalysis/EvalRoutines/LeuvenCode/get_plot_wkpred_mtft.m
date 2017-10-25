function dataout = get_plot_wkpred_mtft(kernels,carrier_freq)
% Specify carrier frequency in Hz. Use carrier_freq = [] for noise carrier


%% Analysis parameters, for cropping/cleaning Wiener kernels
minlag = 0.001;
maxlag = 0.011;
n_eigs = 13;        % Only the first n_eigs eigenvector are retained in h2, to reduce noise 
inspect_eigs = 1;
psth_plots = 1;
fs = 25000;         % Sampling frequency of kernels


%% Stimulus parameters
fm_lo = 4;                              % Minimum modulation frequency in Hz
fm_hi = 1025;
steps_per_octave = 3;
stim_dur = 1;                           % Stimulus duration in sec
ramp_dur = 0.01;                        % Ramps are cosine square
noiseband = [100 10000];                % Frequency band of noise carrier
b_bp = fir1(4000,[noiseband]/fs*2);     % Band pass filter for noise carrier
spl = 50;                               % dB SPL of stimulus


%% Vector of AM frequencies
i = 1;
fm = fm_lo;
while fm < fm_hi
    all_fms(i) = round(fm);
    fm = fm * 2^(1/steps_per_octave); % next freq
    i = i+1;
end
all_fms_4plot = [2 all_fms];
all_fms = [0 all_fms];            % Add 0 Hz modulation frequency
nfms = length(all_fms);
dataout.modfreq = all_fms_4plot;

%% Crop wiener kernels
t = kernels.t;
kernels.h1(t>maxlag) = [];
kernels.h1(t<minlag) = [];
kernels.h2(t>maxlag,:) = [];
kernels.h2(:,t>maxlag) = [];
kernels.h2(t<minlag,:) = [];
kernels.h2(:,t<minlag) = [];


%% Decompose into eigenvectors
[uENV,sENV] = eig(kernels.h2);
sENV = diag(sENV);
for i = 1:24
    eigvec(i,:) = uENV(:,i)'*sENV(i);
end


%% Reconstruct h2 based on first n_eigs eigenvectors. This reduces noise
kernels.h2 = 0*kernels.h2;
for i = 1:(min([n_eigs length(sENV)]))
    kernels.h2 = kernels.h2 + uENV(:,i)*uENV(:,i)'*sENV(i);
end


%% Build carrier signal and gate
stim_pts = round(stim_dur*fs);
t = (0:(stim_pts-1))/fs;
stim_gate = tukeywin(stim_pts,2*ramp_dur/stim_dur)'; %raised cosine ramps
if ~isempty(carrier_freq)                               % Tone MTF
    carrier_sig = cos(2*pi*carrier_freq*t);             % Tone
    desired_rms_Pa = (10^(spl/20))*20e-6;               % Desired RMS amplitude of stimuli in Pa
else                                                    % Noise MTF
    spec_level = spl - 10*log10(diff(noiseband));       % Spectrum level of noise
    noise_spl = spec_level + 10*log10(fs/2);            % SPL of white noise (BW: 0 - fn) with desired spectrum level
    noise_rms = 10^(noise_spl/20)*20e-6;                % RMS amplitude
    randn('state',0);                                   % Frozen noise
    BBN_Pa = noise_rms*randn(stim_pts,1)';              % White noise (BW: 0 - fn) with appropriate spectrum level
    carrier_sig = conv(BBN_Pa,b_bp,'same');             % Band-limited noise of appropriate spectrum level
    desired_rms_Pa = rms(carrier_sig);                  % RMS amplitude of noise waveform prior to modulation (target RMS)
end


%% Main loop
tic
progress_bar = waitbar(0,'Running simulations...');
for i = 1:nfms
    waitbar(i/nfms,progress_bar)

    
    %% Build stimulus
    modulator = 1 + cos(2*pi*all_fms(i)*t);                 % modulator
    stim = carrier_sig.*modulator;                          % AM tone
    stim_Pa = stim * desired_rms_Pa / rms(stim);            % In Pa
    stim_Pa = stim_Pa.*stim_gate;                           % Gated
    
    
    %% Get prediction
    psth = wkpred(stim_Pa,kernels,fs);
    pred_h1(i,:) = (psth.h0+psth.h1);
     % pred_h2(i,:) = (psth.h0+psth.h2);
        pred_h2(i,:) = (psth.h2);
    pred_all(i,:) = (psth.h0+psth.h1+psth.h2);
    
    
    %% Calculate response synchrony to mod freq
    foi = all_fms(i);                                       % frequency of interest
    if foi == 0                                             % ...the unmodulated stimulus
        foi = 100;
    end
    analysis_cycles = (stim_dur - 2*ramp_dur)*foi;          % Ensure integer number of cycles
    analysis_window = round(ramp_dur*fs):round((ramp_dur+analysis_cycles/foi)*fs);
    nfft = fs;                                              % For 1 Hz resolution
    pred_windowed = pred_h2(i,analysis_window);             % Focus is on h2 becuase h1 doesn't show synch to env
    L = length(pred_windowed);
    pred_fft = fft(pred_windowed,nfft)/L;
    pred_amp = 2*abs(pred_fft(1:nfft/2+1));
    pred_amp = pred_amp/pred_amp(1,1)/2;
    f = fs/2*linspace(0,1,nfft/2+1);
    [~,foi_ind] = min(abs(f-foi));
    vs(i) = pred_amp(foi_ind);
    
    
    rate.h1(i) = rms(pred_h1(i,:));
    rate.h2(i) = mean(pred_h2(i,:));
    rate.all(i) = mean(pred_all(i,:));
end
toc
close(progress_bar)

dataout.rate = rate.h2;
dataout.synch = vs;

return

%% Figures~!
if inspect_eigs
    figure
    maxamp = max(max(abs(eigvec)));
    for i = 1:24
        subplot(6,4,i)
        plot(eigvec(i,:))
        ylim([-maxamp*1.1 maxamp*1.1])
    end
end

if psth_plots
    figure
    ncolumns = 3;
    nrows = ceil(nfms/ncolumns);
    ymax = max(max(pred_all));
    for iplot = 1:nfms
        subplot(nrows,ncolumns,iplot)
        plot(t,pred_h1(iplot,:),'r-',t,pred_h2(iplot,:),'b-')
        ylim([0 ymax*1.1])
    end
    
end

figure
% subplot(1,2,1)
semilogx(all_fms_4plot,rate.h2,'.b-')
set(gca,'XLim',[1 max(all_fms)*2],'XTick',[2 4 8 16 32 64 128 256 512 1024],...
    'XTickLabel',[' 0  ';' 4  ';' 8  ';' 16 ';' 32 ';' 64 ';' 128';' 256';' 512';'1024'],...
    'YLim',[0 1.1*max(rate.all)])

figure
semilogx(all_fms_4plot,vs,'.b-')
set(gca,'XLim',[1 max(all_fms)*2],'XTick',[2 4 8 16 32 64 128 256 512 1024],...
    'XTickLabel',[' 0  ';' 4  ';' 8  ';' 16 ';' 32 ';' 64 ';' 128';' 256';' 512';'1024'],...
    'YLim',[0 1])

ddd = 1;
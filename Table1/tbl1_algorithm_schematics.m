%% Configure Library Path
util_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/utils');
addpath(util_path);
%% Set Colors
orange = '#fc6600';
purple = '#a717be';
black = '#031B28';
brown = '#6e2f1a';
green = '#c7d9a5';
blue = '#00568E';
red = '#D0051A';
red_rgb = [208,5,26]/255;
%% Simulate Signals
% [1] Set Parameters
method = 'tukey';
T = 5;
f1 = 25; f2 = 40;
nBurst1 = 10;
nBurst2 = 10;
noise_snr = 5;
rng_pks = [5 6];
rnd_seed = 97;
sim_verbose = false;
% [2] Simulate a Signal
[tvec,simulated_signal,Fs,data,btime_f1,burst_f1] = generate_simulation_rnd(method,T,f1,f2,nBurst1,nBurst2,noise_snr,rng_pks,rnd_seed,sim_verbose);
%% Plot Time Domain Analysis: Bandpass Filtering Method (BP)
% [1] Filter Signal
lf = f1-2; hf = f1+2; 
filt_signal = eegfilt(simulated_signal',Fs,lf,hf);
% [2] Detect Bursts
thrA_factor = 1.2;
[btime_ts, burst_ts] = detect_burst_timeseries(Fs, tvec, filt_signal, lf, hf, thrA_factor);
thr_ts = mean(filt_signal) + thrA_factor*std(filt_signal);
subthr_ts = mean(filt_signal) + (thrA_factor-1)*std(filt_signal);
durthr_ts = 2*(Fs/((lf+hf)/2));
% [3] Plot Filtered Signal
figure();
plot(tvec,filt_signal,'Color',blue,'LineWidth',2);
box off;
set(gca,'TickLength',[0 0],'XColor','none','YColor','none','Color','none');
set(gcf,'Units','Normalized','OuterPosition',[0.13,0.11,0.325,0.425],'Color','w');
xlim([3.5,4.5]);
line([3.9,4.1],[-1.2,-1.2],'LineWidth',2.5,'Color',black);
line([4,4],[-1.2,-0.7],'LineWidth',2.5,'Color',black);
text(3.92,-1.4,'200ms','FontSize',26,'FontName','Helvetica');
text(4.03,-0.95,'0.5mV','FontSize',26,'FontName','Helvetica');
% [4] Plot Detected Bursts
figure(); hold on;
set(gca,'DefaultLineLineWidth',2.0);
set(gcf,'Units','Normalized','OuterPosition',[0.13,0.11,0.325,0.425],'Color','w');
for i = 1:length(btime_ts)
    btime_mid = btime_ts{i}(round(length(btime_ts{i})/2));
    rectangle('Position',[btime_mid-(durthr_ts/(Fs*2)),min(filt_signal),durthr_ts/(Fs),abs(min(filt_signal))+max(filt_signal)],'EdgeColor',green,'FaceColor',green); % plot duration threshold
end
plot(tvec,filt_signal,'Color',blue);
yln_ts = yline(thr_ts,'Color',orange,'LineWidth',3,'LineStyle','-.','alpha',1);
yln_ts.LabelVerticalAlignment = 'top';
yln_ts.LabelHorizontalAlignment = 'right';
yln_ts2 = yline(subthr_ts,'Color',orange,'LineWidth',3,'LineStyle','-','alpha',1);
yln_ts2.LabelVerticalAlignment = 'bottom';
yln_ts2.LabelHorizontalAlignment = 'right';
for i = 1:length(btime_ts)
    plot(btime_ts{i},burst_ts{i},'Color',red);
end
set(gca,'TickLength',[0 0],'XColor','none','YColor','none','Color','none');
set(gcf,'Color','w');
xlim([3.5,4.5]);
% [5] Plot Zoom-In Visualization
figure(); hold on;
set(gcf,'Units','Normalized','OuterPosition',[0.13,0.31,0.275,0.315],'Color','w');
set(gca,'DefaultLineLineWidth',2.0);
start = find(tvec >= btime_ts{9}(1)-0.1,1,'first');
stop = find(tvec <= btime_ts{9}(end)+0.1,1,'last');
plot(tvec(start:stop),filt_signal(start:stop),'Color',blue);
plot(btime_ts{9},burst_ts{9},'Color',red);
xlim([tvec(start)-0.1 tvec(stop)+0.1]);
set(gca,'TickLength',[0 0],'XColor','none','YColor','none','Color','none');
set(gcf,'Color','w'); hold off;
xlim([btime_ts{9}(1)-0.2,btime_ts{9}(end)+0.2]);
line([tvec(stop)-0.15,tvec(stop)-0.05],[-1,-1],'LineWidth',2.5,'Color',black);
line([tvec(stop)-0.05,tvec(stop)-0.05],[-1,-0.5],'LineWidth',2.5,'Color',black);
text(tvec(stop)-0.15,-1.15,'100ms','FontSize',23,'FontName','Helvetica');
text(tvec(stop)-0.04,-0.75,'0.5mV','FontSize',23,'FontName','Helvetica');
% [6] Plot Manual Legend Box
leg_fig = figure(); hold on;
set(gcf,'Color','w');
set(gca,'DefaultLineLineWidth',10.0,'Color','none');
plot(rand(1,2),rand(1,2),'Color',blue);
plot(rand(1,2),rand(1,2),'Color',red);
plot(rand(1,2),rand(1,2),'Color',orange,'LineStyle',':');
plot(rand(1,2),rand(1,2),'Color',orange);
plot(rand(1,2),rand(1,2),'Color',green);
lgnd = legend('Filtered Signal','Detected Burst','Amplitude Threshold','Amplitude Sub-threshold','Duration Threshold','FontSize',34,'Color','w','NumColumns',1);
lgnd.LineWidth = 4;
lgnd.ItemTokenSize = [50,8];
makeLegendToImage(leg_fig,lgnd,'line');
%% Plot Time Domain Analysis: Envelope-Based Method (ENV)
% [1] Detect Bursts
env_signal = abs(hilbert(filt_signal));
prcthrA = 70;
[btime_env,burst_env,benvelope] = detect_burst_ampenv(Fs,tvec,filt_signal,lf,hf,prcthrA);
thr_env = prctile(env_signal,prcthrA);
durthr_env = 2*(Fs/((lf+hf)/2));
% [2] Plot Filtered Signal and Envelope
figure(); hold on;
set(gcf,'Units','Normalized','OuterPosition',[0.13,0.11,0.325,0.425],'Color','white');
set(gca,'TickLength',[0 0],'XColor','none','YColor','none','Color','none');
plot(tvec,filt_signal,'Color',blue,'LineWidth',2.0);
plot(tvec,env_signal,'Color',purple,'LineWidth',3.0);
box off;
xlim([3.5,4.5]);
line([3.9,4.1],[-1.2,-1.2],'LineWidth',2.5,'Color',black);
line([4,4],[-1.2,-0.7],'LineWidth',2.5,'Color',black);
text(3.92,-1.4,'200ms','FontSize',26,'FontName','Helvetica')
text(4.03,-0.95,'0.5mV','FontSize',26,'FontName','Helvetica')
% [3] Plot Detected Bursts
figure(); hold on;
set(gcf,'Units','Normalized','OuterPosition',[0.13,0.11,0.325,0.425],'Color','white');
env_limited = env_signal;
env_limited(env_limited >= thr_env) = thr_env;
fill(tvec,env_limited,'white','LineStyle','none');
for i = 1:length(btime_env)
    btime_mid = btime_env{i}(round(length(btime_env{i})/2));
    rectangle('Position',[btime_mid-(durthr_env/(Fs*2)),min(env_signal),durthr_env/(Fs),abs(min(env_signal))+max(env_signal)],'EdgeColor',green,'FaceColor',green); % plot duration threshold
    benvelope{i}(1) = thr_env; % fill the gap on the front
    benvelope{i}(end) = thr_env; % fill the gap on the end
    fill(btime_env{i},benvelope{i},red_rgb,'LineStyle','none','FaceAlpha',0.45);
end
plot(tvec,env_signal,'Color',purple,'LineWidth',3.0);
yln_env = yline(thr_env,'Color',orange,'LineWidth',3,'LineStyle','--','alpha',1);
yln_env.LabelVerticalAlignment = 'top';
yln_env.LabelHorizontalAlignment = 'right';
set(gca,'TickLength',[0 0],'XColor','none','YColor','none','Color','none');
xlim([3.5,4.5]);
% [4] Plot Zoom-In Visualization
start = find(tvec >= btime_env{6}(1)-0.1,1,'first');
stop = find(tvec <= btime_env{6}(end)+0.1,1,'last');
figure(); hold on;
set(gcf,'Units','Normalized','OuterPosition',[0.13,0.31,0.257,0.345],'Color','white');
set(gca,'DefaultLineLineWidth',1.5,'TickLength',[0 0],'XColor','none','YColor','none','Color','none');
fill(btime_env{6},benvelope{6},red_rgb,'LineStyle','none','FaceAlpha',0.6);
fill(tvec(start:stop),env_limited(start:stop),'white','LineStyle','none');
plot(tvec(start:stop),env_signal(start:stop),'Color',purple,'LineWidth',2);
xlim([btime_env{6}(1)-0.2,btime_env{6}(end)+0.2]);
line([tvec(stop)-0.1,tvec(stop)],[-0.2,-0.2],'LineWidth',2.5,'Color',black);
line([tvec(stop),tvec(stop)],[-0.2,-0],'LineWidth',2.5,'Color',black);
text(tvec(stop)-0.1,-0.26,'100ms','FontSize',26,'FontName','Helvetica');
text(tvec(stop)+0.01,-0.1,'0.2mV','FontSize',26,'FontName','Helvetica');
% [5] Plot Manual Legend Box
leg_fig = figure(); hold on;
set(gcf,'Color','w');
set(gca,'DefaultLineLineWidth',10.0,'Color','none');
p1 = plot(rand(1,2),rand(1,2),'Color',blue);
p2 = plot(rand(1,2),rand(1,2),'Color',purple);
p3 = plot(rand(1,2),rand(1,2),'Color',red);
p4 = plot(rand(1,2),rand(1,2),'Color',orange,'LineStyle',':');
p5 = plot(rand(1,2),rand(1,2),'Color',green);
lgnd = legend('Filtered Signal','Envelope Signal','Detected Burst','Amplitude Threshold','Duration Threshold','FontSize',34,'Color','w','NumColumns',1);
lgnd.LineWidth = 4;
lgnd.ItemTokenSize = [50,8];
makeLegendToImage(leg_fig,lgnd,'line');
%% Set Visualization Parameters for Spectral Analysis-Based Algorithms
fnt_size = 68;
outer_pos = [1.023,-0.047,1.108,1.722];
inner_pos = [0.13,0.20,0.433,0.632];
inner_pos_extend = [0.13,0.20,0.433,0.632];
%% Compute Spectrograms
% [1] Single-Tapered Short Time Fourier Transform (S-STFT)
freq_range = unique([15:5:50 50:10:70]);
dt = 0.01; nOscCycle = 6;
verbose = false; interp_opt = false;
[Spec_f_stp,Spec_t_stp,Spec_stp,fft_interval_stp] = spectrogram_stp(tvec, simulated_signal, Fs, freq_range, dt, nOscCycle, verbose, interp_opt);
% [2] Multitaper Spectrogram (MTP)
nw = 2; ntp = 3;
plot_psd = false; check_opt = true; plot_tp = false;
[Spec_f_mtp,Spec_t_mtp,Spec_mtp,fft_interval_mtp] = spectrogram_mtp(tvec, simulated_signal, Fs, nw, ntp, freq_range, dt, nOscCycle, interp_opt, ... 
    plot_psd, check_opt, plot_tp, verbose);
% [3] Continuous Wavelet Transform (CWT)
wname = 'amor'; scale_opt = 'default';
[wav_freq,wav_time,Spec_cwt,coi] = spectrogram_cwt(tvec, simulated_signal, Fs, wname, scale_opt, verbose);
%% Plot Time-Frequency Domain Analysis: S-STFT
% [1] Detect Bursts
thrA_factor = 1.0;
[stime_stp,sburst_stp,binSpec_stp] = detect_burst_spectrogram(Spec_f_stp,Spec_t_stp,Spec_stp,dt,thrA_factor,verbose);
% [2] Interpolation
interp_method = 'spline';
Spec_f_ip = Spec_f_stp(1):Spec_f_stp(end);
Spec_t_ip = Spec_t_stp;
[Xq, Yq] = meshgrid(Spec_f_ip,Spec_t_ip);
% visualized center frequency may deviate from the actual center
% frequency after interpolation
SpecInterp = interp2(Spec_f_stp,Spec_t_stp,Spec_stp,Xq,Yq,interp_method);
% [3] Plot Power Spectrogram
figure(); box off;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(Spec_t_ip,Spec_f_ip,imgaussfilt(SpecInterp',2)); % spectral power may decrease when Gaussian filtered
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
c = colorbar('Location','northoutside');
ylabel(c, 'Power (mV^2)','FontSize',fnt_size);
xlim([3.5 4.5]); ylim([15 45]);
set(gca,'TickLength',[0 0],'XColor','none','InnerPosition',inner_pos_extend,'LineWidth',0.1,'FontSize',fnt_size);
hold on;
line([4.19 4.49],[15.8,15.8],'Color','w','LineWidth',20);
text(4.21,17.7,'300ms','FontSize',60,'FontName','Helvetica','Color','w');
% [4] Plot Quantized Spectrogram
figure(); box on;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(Spec_t_stp,Spec_f_stp,imcomplement(binSpec_stp'));
axis xy; colormap(gray);
xlim([3.5 4.5]); ylim([15 45]);
ylabel('Frequency (Hz)');
line([stime_stp{3,3}(1) stime_stp{3,3}(end)], [Spec_f_stp(3) Spec_f_stp(3)],'Color',red,'LineWidth',20,'LineStyle','-');
line([stime_stp{3,3}(1) stime_stp{3,3}(1)], [Spec_f_stp(3)-1 Spec_f_stp(3)+1],'Color',red,'LineWidth',20,'LineStyle','-');
line([stime_stp{3,3}(end) stime_stp{3,3}(end)], [Spec_f_stp(3)-1 Spec_f_stp(3)+1],'Color',red,'LineWidth',20,'LineStyle','-');
set(gca,'TickLength',[0 0],'XTickLabel',{[]},'InnerPosition',inner_pos,'LineWidth',3,'FontSize',fnt_size);
% [5] Plot Zoom-In Visualization
figure(); box off;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(Spec_t_ip,Spec_f_ip,imgaussfilt(SpecInterp',2)); % spectral power may decrease when Gaussian filtered
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%c = colorbar('Location','northoutside'); ylabel(c, 'Power (mV^2/Hz)','FontSize',fnt_size);
xlim([3.6 3.9]); ylim([17 32]);
line([stime_stp{3,3}(1) stime_stp{3,3}(end)], [fft_interval_stp(3,1) fft_interval_stp(3,1)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_stp{3,3}(1) stime_stp{3,3}(end)], [fft_interval_stp(3,2) fft_interval_stp(3,2)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_stp{3,3}(1) stime_stp{3,3}(1)], [fft_interval_stp(3,1) fft_interval_stp(3,2)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_stp{3,3}(end) stime_stp{3,3}(end)], [fft_interval_stp(3,1) fft_interval_stp(3,2)],'Color',red,'LineWidth',10,'LineStyle','--');
set(gca,'TickLength',[0 0],'XColor','none','InnerPosition',inner_pos,'LineWidth',1,'FontSize',fnt_size);
%text(stime_stp{3,2}(1)-0.06,fft_interval_stp(3,2)+1,'Above Threshold','Color',red,'FontSize',fnt_size);
line([3.79 3.89],[17.5,17.5],'Color','w','LineWidth',20);
text(3.8,18.6,'100ms','FontSize',60,'FontName','Helvetica','Color','w');
% [6] Plot Manual Legend Box
leg_fig = figure(); hold on;
set(gcf,'Color','w');
set(gca,'DefaultLineLineWidth',10.0,'Color','none');
plot(rand(1,2),rand(1,2),'Color',red);
plot(rand(1,2),rand(1,2),'Color',black);
lgnd = legend({'Detected Burst',['Points Above' newline 'Amplitude Threshold']},'FontSize',34,'Color','w','NumColumns',1);
lgnd.LineWidth = 4;
lgnd.ItemTokenSize = [50,8];
makeLegendToImage(leg_fig,lgnd,'line');
%% Plot Time-Frequency Domain Analysis: MTP
% [1] Detect Bursts
[stime_mtp,sburst_mtp,binSpec_mtp] = detect_burst_spectrogram(Spec_f_mtp,Spec_t_mtp,Spec_mtp,dt,thrA_factor,verbose);
% [2] Interpolation
interp_method = 'spline';
Spec_f_ip = Spec_f_mtp(1):Spec_f_mtp(end);
Spec_t_ip = Spec_t_mtp;
[Xq, Yq] = meshgrid(Spec_f_ip,Spec_t_ip);
SpecInterp = interp2(Spec_f_mtp,Spec_t_mtp,Spec_mtp,Xq,Yq,interp_method);
% [3] Plot Power Spectrogram
figure(); box off;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(Spec_t_ip,Spec_f_ip,imgaussfilt(SpecInterp',2));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
c = colorbar('Location','northoutside');
ylabel(c, 'Power (mV^2)','FontSize',fnt_size);
xlim([1.9 2.9]); ylim([15 45]);
set(gca,'TickLength',[0 0],'XColor','none','InnerPosition',inner_pos_extend,'LineWidth',0.1,'FontSize',fnt_size);
hold on;
line([2.59 2.89],[15.8,15.8],'Color','w','LineWidth',20);
text(2.61,17.7,'300ms','FontSize',60,'FontName','Helvetica','Color','w');
% [4] Plot Quantized Spectrogram
figure(); box on;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(Spec_t_mtp,Spec_f_mtp,imcomplement(binSpec_mtp'));
axis xy; colormap(gray);
xlim([1.9 2.9]); ylim([15 45]);
ylabel('Frequency (Hz)');
line([stime_mtp{3,2}(1) stime_mtp{3,2}(end)], [Spec_f_mtp(3) Spec_f_mtp(3)],'Color',red,'LineWidth',20,'LineStyle','-');
line([stime_mtp{3,2}(1) stime_mtp{3,2}(1)], [Spec_f_mtp(3)-1 Spec_f_mtp(3)+1],'Color',red,'LineWidth',20,'LineStyle','-');
line([stime_mtp{3,2}(end) stime_mtp{3,2}(end)], [Spec_f_mtp(3)-1 Spec_f_mtp(3)+1],'Color',red,'LineWidth',20,'LineStyle','-');
set(gca,'TickLength',[0 0],'XTickLabel',{[]},'InnerPosition',inner_pos,'LineWidth',3,'FontSize',fnt_size);
% [5] Plot Zoom-In Visualization
figure(); box off;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(Spec_t_ip,Spec_f_ip,imgaussfilt(SpecInterp',2)); % spectral power may decrease when Gaussian filtered
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
xlim([2.2 2.7]); ylim([17 32]);
line([stime_mtp{3,2}(1) stime_mtp{3,2}(end)], [fft_interval_mtp(3,1) fft_interval_mtp(3,1)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_mtp{3,2}(1) stime_mtp{3,2}(end)], [fft_interval_mtp(3,2) fft_interval_mtp(3,2)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_mtp{3,2}(1) stime_mtp{3,2}(1)], [fft_interval_mtp(3,1) fft_interval_mtp(3,2)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_mtp{3,2}(end) stime_mtp{3,2}(end)], [fft_interval_mtp(3,1) fft_interval_mtp(3,2)],'Color',red,'LineWidth',10,'LineStyle','--');
set(gca,'TickLength',[0 0],'XColor','none','InnerPosition',inner_pos,'LineWidth',1,'FontSize',fnt_size);
line([2.58 2.68],[17.5,17.5],'Color','w','LineWidth',20);
text(2.56,18.6,'100ms','FontSize',60,'FontName','Helvetica','Color','w');
%% Plot Time-Frequency Domain Analysis: Continuous Wavelet Transform (CWT)
% [1] Detect Bursts
dt_cwt = 1/Fs;
[stime_cwt,sburst_cwt,binSpec_cwt] = detect_burst_spectrogram(wav_freq,wav_time,Spec_cwt',dt_cwt,thrA_factor,verbose);
fidx = find(abs(wav_freq-f1) == min(abs(wav_freq-f1)));
foi = wav_freq(fidx);
interval_cwt = [foi-(foi/(6*sqrt(2)))/2 foi+(foi/(6*sqrt(2)))/2];
% [2] Crop Wavelet Spectrogram
crop_start = find(wav_freq < 50,1,'first');
crop_end = find(wav_freq > 15,1,'last');
wav_freq = wav_freq(crop_start:crop_end);
Spec_cwt = Spec_cwt(crop_start:crop_end,:);
binSpec_cwt = binSpec_cwt(:,crop_start:crop_end);
% [3] Plot Power Spectrogram
figure(); hold on; box off;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(wav_time,wav_freq,imgaussfilt(Spec_cwt,1));
shading flat;
plot(wav_time,coi,'--w','LineWidth',2);
xlim([3 4]);
ylim([wav_freq(end) wav_freq(1)]);
ftemp = flip(wav_freq);
ax = gca;
ax.YTick = round(ftemp(1:5:end));
set(ax,'yscale','log');
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
c = colorbar('Location','northoutside');
ylabel(c,'Coefficient Magnitude','FontSize',fnt_size);
set(gca,'TickLength',[0 0],'XColor','none','InnerPosition',inner_pos_extend,'LineWidth',0.1,'FontSize',fnt_size);
hold on;
line([3.69 3.99],[16.1,16.1],'Color','w','LineWidth',20);
text(3.71,17.5,'300ms','FontSize',60,'FontName','Helvetica','Color','w');
% [4] Plot Quantized Spectrogram
figure(); box on;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(wav_time,wav_freq,imcomplement(binSpec_cwt'));
ftemp = flip(wav_freq);
ax = gca;
ax.YTick = round(ftemp(1:5:end));
set(ax,'yscale','log');
axis xy; colormap(gray);
xlim([3 4]); ylim([15 45]);
ylabel('Frequency (Hz)');
line([stime_cwt{fidx,3}(1) stime_cwt{fidx,3}(end)], [foi foi],'Color',red,'LineWidth',20,'LineStyle','-');
line([stime_cwt{fidx,3}(1) stime_cwt{fidx,3}(1)], [foi-1 foi+1],'Color',red,'LineWidth',20,'LineStyle','-');
line([stime_cwt{fidx,3}(end) stime_cwt{fidx,3}(end)], [foi-1 foi+1],'Color',red,'LineWidth',20,'LineStyle','-');
set(gca,'TickLength',[0 0],'XTickLabel',{[]},'InnerPosition',inner_pos,'LineWidth',3,'FontSize',fnt_size);
% [5] Plot Zoom-In Visualization
figure(); hold on; box off;
set(gcf,'Units','Normalized','OuterPosition',outer_pos,'Color','w');
imagesc(wav_time,wav_freq,imgaussfilt(Spec_cwt,1));
shading flat;
plot(wav_time,coi,'--w','LineWidth',2);
xlim([3.58 3.92]); ylim([19 37]);
ftemp = flip(wav_freq);
ax = gca;
ax.YTick = round(ftemp(1:3:end));
set(ax,'yscale','log');
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'TickLength',[0 0],'XColor','none','InnerPosition',inner_pos_extend,'LineWidth',0.1,'FontSize',fnt_size);
line([stime_cwt{fidx,3}(1) stime_cwt{fidx,3}(end)], [interval_cwt(1) interval_cwt(1)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_cwt{fidx,3}(1) stime_cwt{fidx,3}(end)], [interval_cwt(2) interval_cwt(2)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_cwt{fidx,3}(1) stime_cwt{fidx,3}(1)], [interval_cwt(1) interval_cwt(2)],'Color',red,'LineWidth',10,'LineStyle','--');
line([stime_cwt{fidx,3}(end) stime_cwt{fidx,3}(end)], [interval_cwt(1) interval_cwt(2)],'Color',red,'LineWidth',10,'LineStyle','--');
line([3.81 3.91],[19.5,19.5],'Color','w','LineWidth',20);
text(3.815,20.6,'100ms','FontSize',60,'FontName','Helvetica','Color','w');
%% Configure Library Path
util_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/data/simulation_data');
addpath(util_path);
addpath(data_path);
%% Set Colors
red = '#d13b3b';
orange = '#d98100';
bluel = '#377A98';
blued = '#244A80';
black = '#242621';
%% Simulate Example Signal
% [1] Set Parameters
method = 'tukey';
T = 10;
f1 = 25; f2 = 40;
nBurst1 = 10;
nBurst2 = 15;
noise_snr = 3;
rng_pks = [3 12];
rnd_seed = 3;
sim_verbose = false;
% [2] Simulate a Signal
[tvec,simulated_signal,Fs,data,btime_f1,burst_f1,cyc1] = generate_simulation_rnd(method,T,f1,f2,nBurst1,nBurst2,noise_snr,rng_pks,rnd_seed,sim_verbose);
%% Plot Simulated Signal Schematics
fnt_sz = 54;
draw_schematics(tvec,data(:,2),{'β Bursts';'(25Hz)'},fnt_sz,bluel);
draw_schematics(tvec,data(:,3),{'γ Bursts';'(40Hz)'},fnt_sz,blued);
draw_schematics(tvec,data(:,4),{'Aperiodic';'Noise'},fnt_sz,orange);
draw_schematics(tvec,data(:,5),{'Gaussian';'Noise'},fnt_sz,red);
draw_schematics(tvec,simulated_signal,{'Simulated';'Signal'},fnt_sz,black,true);
%% Plot Welch's Power Spectral Density Estimate (Entire Dataset)
method = 'tukey';
T = 300;           % unit: s
f1 = 25; f2 = 40;  % unit: Hz
nBurst1 = 30;
nBurst2 = 30;
rng_pks = [3 12];
sim_verbose = false;
% [1] Get Random Parameters
nSimulationPerNoise = 1000; % number of simulations per each SNR
listNoise = -10:2:10;
if isfile('randseed.mat')
    listSeed = load('randseed.mat').randseed;
else
    listSeed = zeros(length(listNoise),nSimulationPerNoise);
    for n = 1:length(listNoise)
        listSeed(n,:) = randsample(1:2^15,nSimulationPerNoise);
    end
end
nSeed = numel(listSeed);
% [2] Simulate Signals
[psd_sig_data, psd_noise_data] = deal(cell(size(listSeed)));
for n = 1:length(listNoise)
    SNR_wn = listNoise(n);
    for s = 1:size(listSeed,2)
        [tvec,simulated_signal,Fs,data,~,~,~] = generate_simulation_rnd(method,T,f1,f2,nBurst1,nBurst2,noise_snr,rng_pks,listSeed(s),sim_verbose);
        psd_window = Fs*30;
        noverlap = floor(psd_window/2); % 50 percent overlap
        nfft = psd_window;
        [pxx_sig,f_sig] = pwelch(simulated_signal,hanning(psd_window),noverlap,nfft,Fs);
        [pxx_noise,f_noise] = pwelch(data(:,4)+data(:,5),hanning(psd_window),noverlap,nfft,Fs);
        psd_sig_data{n,s} = pxx_sig;
        psd_noise_data{n,s} = pxx_noise;
    end
end
psd_sig = sum([psd_sig_data{:}],2)./nSeed;
psd_noise = sum([psd_noise_data{:}],2)./nSeed;
% [3] Plot Welch's PSD Estimate
figure(); hold on;
plot(f_sig,psd_sig,'Linewidth',5,'Color','#242621');
plot(f_noise,psd_noise,'Linewidth',5,'Color','#d98100');
xlim([-0.5 60]);
xlabel('Frequency (Hz)');
ylabel('PSD (mV^2/Hz)');
xticks(0:10:60);
yticks(0:0.005:0.03);
ylim([0,0.026]);
lgnd = legend('Simulated Signal','Simulated Noise');
lgnd.FontSize = 42;
set(gca,'LineWidth',6,'TickDir','out','FontSize',54,'FontName','Helvetica','FontWeight','bold');
set(gcf,'Units','Normalized','Color','w','Position',[0.205729166666667,0,0.6546875,0.9037]);
%% Plot Welch's Power Spectral Density Estimate (Different SNRs)
% [1] Plot Welch's PSD Estimate
psd_sig_snr = cell(1,length(listNoise));
seed = 1997;
for n = 1:length(listNoise)
    [tvec,simulated_signal,Fs,data,~,~,~] = generate_simulation_rnd(method,T,f1,f2,nBurst1,nBurst2,listNoise(n),rng_pks,seed,sim_verbose);
    psd_window = Fs*10;
    noverlap = floor(psd_window/2); % 50 percent overlap
    nfft = psd_window;
    [pxx_sig,f_sig] = pwelch(simulated_signal,hanning(psd_window),noverlap,nfft,Fs);
    psd_sig_snr{n} = pxx_sig;
end
figure(); hold on;
xl{1} = xline(f1,'--',[num2str(f1) ' Hz'],'Color',bluel);
xl{2} = xline(f2,'--',[num2str(f2) ' Hz'],'Color',blued);
colororder(hot(length(listNoise)+5));
[p, leg_lbl] = deal(cell(1,length(listNoise)));
for n = 1:length(listNoise)
    p{n} = plot(f_sig,psd_sig_snr{n},'LineWidth',3);
    leg_lbl{n} = [num2str(listNoise(n)) ' dB'];
end
xlim([-0.5 60]);
xlabel('Frequency (Hz)');
ylabel('PSD (mV^2/Hz)');
xticks(0:10:60);
yticks(0:0.01:0.07);
ylim([0,0.05]);
% (A) Set Legend
leg = legend([p{:}],leg_lbl);
title(leg,'SNR Levels');
% (B) Figure Settings
set([xl{:}],'LineWidth',3.5,'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','FontSize',48);
set(leg,'FontSize',26,'Position',[0.8359,0.1866,0.1193,0.4671]);
set(gca,'LineWidth',6,'TickDir','out','FontName','Helvetica','FontSize',54,'FontWeight','bold','Position',[0.1846,0.1875,0.6014,0.7500]);
set(gcf,'Units','Normalized','Color','w','Position',[0.205729166666667,0,0.6546875,0.9037]);
%% Appendix: In-Script Functions
% Function #1: Plot Simulation Schematics
function draw_schematics(tvec,signal,label,fnt_sz,clr,scale_opt)
    if nargin < 6
        scale_opt = false;
    end
    figure(); hold on;
    plot(tvec,signal,'Color',clr,'LineWidth',4);
    text(6.52,0,label,'FontSize',fnt_sz,'FontWeight','bold','HorizontalAlignment','center');
    ylim([-3.0,2.3]);
    gcf_position = [0.1328125,0.519444444444444,0.706770833333333,0.343518518518519];
    if scale_opt
        line([8.8,9.3],[-2.8,-2.8],'LineWidth',6.0,'Color',clr);
        line([9.3,9.3],[-2.8,-1.8],'LineWidth',6.0,'Color',clr);
        text(9.32,-2.3,'1mV','FontSize',fnt_sz,'FontName','Helvetica','FontWeight','bold');
        text(8.86,-3.2,'500ms','FontSize',fnt_sz,'FontName','Helvetica','FontWeight','bold');
    end
    set(gca,'TickLength',[0 0],'XColor','none','YColor','none','xlim',[6.9 9.3],'Color','none', ... 
        'Position',[0.204668075049564,0.11,0.696652286312493,0.82]);
    set(gcf,'Units','Normalized','Position',gcf_position,'Color','white');
end